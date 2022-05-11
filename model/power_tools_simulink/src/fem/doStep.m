function [Tr, n_fem_bodies, n_fem_grids, n_th_bodies, n_fv_grids, n_fv_bodies, fem_body_names, ...
    fem_grid_names, th_body_names, fv_grid_names, fv_body_names, fv_gradient] = doStep(solver_type, ...
    dt, max_sub_iter, T_prev, total_time, max_res, model, number_global_dofs, n_fv_gradient_subiters, fv_gradient, ...
    ele_source_term, bus_input)

coder.allowpcode('plain')
coder.extrinsic('licensing.internal.checkOnlineLicense')
coder.extrinsic('warning');

fem_body_names = fieldnames(model.fem_bodies);
n_fem_bodies = numel(fem_body_names);

fem_grid_names = fieldnames(model.fem_grids);
n_fem_grids = numel(fem_grid_names);

free_dofs = model.free_dofs;
lock_dofs = model.lock_dofs;

if isfield(model, "th_bodies")
    th_body_names = fieldnames(model.th_bodies);
    n_th_bodies = numel(th_body_names);
else
    n_th_bodies = 0;
    th_body_names = [];
end

if isfield(model, "fv_bodies")
    fv_body_names = fieldnames(model.fv_bodies);
    n_fv_bodies = numel(fv_body_names);
else
    n_fv_bodies = 0;
    fv_body_names = [];
end

if isfield(model, "fv_grids")
    fv_grid_names = fieldnames(model.fv_grids);
    n_fv_grids = numel(fv_grid_names);
else
    fv_grid_names = [];
    n_fv_grids = 0;
end

if isfield(model, "contacts_fem_fvm")
    contacts_fem_fvm_names = fieldnames(model.contacts_fem_fvm);
    n_contacts_femfvm = numel(contacts_fem_fvm_names);
else
    n_contacts_femfvm = 0;
end

if isfield(model, "contacts_fvm_fvm")
    contacts_fvm_fvm_names = fieldnames(model.contacts_fvm_fvm);
    n_contacts_fvmfvm = numel(contacts_fvm_fvm_names);
else
    n_contacts_fvmfvm = 0;
end

if isfield(model, "signals")
    signal_names = fieldnames(model.signals);
    n_signals = numel(signal_names);
else
    n_signals = 0;
end

first_simulink_cycle = false;
if total_time == 0
    % check license using licensing methods
    [is_valid, lic_type, secret] = licensing.internal.checkOnlineLicense('d8N3!fqFN6Dp2meW', 'Q-Bat');
    
    if ~(is_valid && strcmp(secret, '#3^Lw#B8p$&Kb6e8'))
        error('License is missing or has expired.');
    end
    
    if ~any(strcmp(lic_type, {'TAH', 'Commercial', 'Commercial Run', 'Research', 'Trial'}))
        error('License is missing or has expired.');
    end
    
    if solver_type ~= 1 && n_fv_bodies > 0
        warning('Cannot use Schur-PCG solver with air. Using Direct solver instead.')
    end
    
    first_simulink_cycle = true;
end

Tr = T_prev;
Tlastr = T_prev;

if ~first_simulink_cycle
    M = zeros(number_global_dofs, number_global_dofs);
    K = zeros(number_global_dofs, number_global_dofs);
    F = zeros(number_global_dofs, 1);
    % nonlinear subiterations
    coder.unroll(false)
    for k = 1 : max_sub_iter
        % add static matrices (fem bodies + fem-fem contacts)
        M = full(sparse(model.M_stat(1, :), model.M_stat(2, :), model.M_stat(3, :), number_global_dofs, number_global_dofs));
        K = full(sparse(model.K_stat(1, :), model.K_stat(2, :), model.K_stat(3, :), number_global_dofs, number_global_dofs));
        F = model.F_stat;
        
        % handle th bodies
        for i=1:n_th_bodies
            lrgr_T = model.th_bodies.(th_body_names{i}).lrgr_T;
            [M_loc, K_loc] = assembleTHBody(model.th_bodies.(th_body_names{i}), T_prev);

            M(lrgr_T,lrgr_T) = M(lrgr_T,lrgr_T) + M_loc;
            K(lrgr_T,lrgr_T) = K(lrgr_T,lrgr_T) + K_loc;
        end
        
        % handle contacts fem-fvm
        for i = 1:n_contacts_femfvm
            [P11, P12, P21, P22] = assembleContactFemFvm(model.contacts_fem_fvm.(contacts_fem_fvm_names{i}), T_prev);
            
            th_body_dofs = model.contacts_fem_fvm.(contacts_fem_fvm_names{i}).th_dofs;
            % th_body_massflow_dofs = model.contacts_fem_fvm.(contacts_fem_fvm_names{i}).th_massflow_dofs;
            fem_body_dofs = model.contacts_fem_fvm.(contacts_fem_fvm_names{i}).fem_dofs;
            
            K(fem_body_dofs, fem_body_dofs) = K(fem_body_dofs, fem_body_dofs) + P11;
            K(fem_body_dofs, th_body_dofs) = K(fem_body_dofs, th_body_dofs) - P12;
            K(th_body_dofs, fem_body_dofs) = K(th_body_dofs, fem_body_dofs) - P21;
            K(th_body_dofs, th_body_dofs) = K(th_body_dofs, th_body_dofs) + P22;
            
        end % for i = 1 : n_contacts_fem_fvm
        
        % contacts fvm-fvm
        for i = 1:n_contacts_fvmfvm
            local_contact = model.contacts_fvm_fvm.(contacts_fvm_fvm_names{i});
            % pass temperature from previous pipe to next pipe
            for m = 1 : n_th_bodies
                if m == local_contact.prev_id
                    prev_dofs_T = model.th_bodies.(th_body_names{m}).lrgr_T;
                    prev_dofs_massflow = model.th_bodies.(th_body_names{m}).lrgr_massflow;
                    for j = 1 : numel(local_contact.next_id)
                        for n = 1 : n_th_bodies
                            if n == local_contact.next_id(j)
                                next_dofs_T = model.th_bodies.(th_body_names{n}).lrgr_T;

                                K(next_dofs_T(1), next_dofs_T(1)) =  K(next_dofs_T(1), next_dofs_T(1)) + local_contact.flow_ratio(j);
                                K(next_dofs_T(1), prev_dofs_T(end)) = K(next_dofs_T(1), prev_dofs_T(1)) - local_contact.flow_ratio(j);

                                % handle mass flows
                                next_dofs_massflow = model.th_bodies.(th_body_names{n}).lrgr_massflow;

                                K(next_dofs_massflow(1), next_dofs_massflow(1)) = 1;
                                K(next_dofs_massflow(1), prev_dofs_massflow(1)) = -local_contact.flow_ratio(j);
                            end
                        end
                    end
                end
            end
        end
        
        % handle bcs
        F_signal = zeros(size(F,1), size(F,2));
        
        Trpr = Tr;
        
        for i=1:n_signals
            local_signal = model.signals.(signal_names{i});
            
            lrgr = local_signal.lrgr;
            
            if ~isempty(local_signal.data_energy)
                temp_data = local_signal.data_energy;
                signal_value = nearestInterp1(temp_data(:,1), temp_data(:,2), bus_input(i));
            else
                signal_value = bus_input(i);
            end
            
            if local_signal.order == 0
                Tr(lrgr(local_signal.localDof)) = signal_value;
                Trpr(lrgr(local_signal.localDof)) = signal_value;
            elseif local_signal.order == 1
                
                if numel(local_signal.Fr) == numel(lrgr)
                    F_signal(lrgr) = F_signal(lrgr) + local_signal.Fr * signal_value;
                end
                
            end
            
        end
        
        % get electric source term
        F_ele = assembleElectroSourceRHS(model.fem_grids, model.fem_bodies, number_global_dofs, ele_source_term);
        
        % get nonlinear parts of fem bodies
        for body_index = 1 : n_fem_bodies
    
            body_name = fem_body_names{body_index};
            local_body = model.fem_bodies.(body_name);
            local_dofs = local_body.lrgr;
            T_body_reduced = Trpr(local_dofs);
            
            g_id = local_body.grid_name;
            
            for j = 1 : n_fem_grids
                if(j == g_id)
                    [M_local, K_local] = getNonstaticMatrices(model.fem_grids.(fem_grid_names{j}), local_body, T_body_reduced);

                    M(local_dofs,local_dofs) = M(local_dofs, local_dofs) + M_local;
                    K(local_dofs,local_dofs) = K(local_dofs, local_dofs) + K_local;
                end
            end
            
        end
        
        % calculate fv correction
        for body_index = 1 : n_fv_bodies
    
            body_name = fv_body_names{body_index};
            local_body = model.fv_bodies.(body_name);
            local_dofs = local_body.lrgr;
            T_body_reduced = Trpr(local_dofs);
    
            grid_name = fv_grid_names{local_body.grid_name};
    
            for j = 1 : n_fv_grids
                if fv_grid_names{j} == grid_name
                    [fv_gradient{body_index}, fv_gradient_faces] = calculateUnsFVGradients(model.fv_grids.(fv_grid_names{j}), T_body_reduced', ...
                        n_fv_gradient_subiters, fv_gradient{body_index});

                    F(local_dofs) = F(local_dofs) + fvsolver.assemblation.assembleUnsFVDiffusionRHS(model.fv_grids.(fv_grid_names{j}).fvmesh, ...
                        model.fv_grids.(fv_grid_names{j}).T, model.fv_grids.(fv_grid_names{j}).lambda, fv_gradient_faces);
                end
            end
        end
        
        F = F + F_ele + F_signal;
        
        A_tmp = M(model.free_dofs, model.free_dofs) / dt + K(model.free_dofs, model.free_dofs);
        F_tmp = F(free_dofs) ...
            + M(free_dofs, free_dofs) / dt * Tlastr(free_dofs) ...
            - 1/dt*M(free_dofs, lock_dofs) * (Trpr(lock_dofs) - Tlastr(lock_dofs))...
            - K(free_dofs, lock_dofs) * Trpr(lock_dofs);
        
        if(solver_type == 1 || n_fv_bodies > 0) % Direct solver
            
            Tr(free_dofs) = A_tmp \ F_tmp;
            
        elseif(solver_type == 2) % Schur-PCG solver
            
            fem_dofs = model.fem_dofs;
            fv_dofs = model.fv_dofs;
            
            A_sch = A_tmp(fem_dofs, fem_dofs);
            B_sch = A_tmp(fem_dofs, fv_dofs);
            C_sch = A_tmp(fv_dofs, fem_dofs);
            D_sch = A_tmp(fv_dofs, fv_dofs);
            
            Fa_sch = F_tmp(fem_dofs);
            Fb_sch = F_tmp(fv_dofs);
            
            z = Fa_sch - B_sch * (D_sch \ Fb_sch);
            
            tmp = (A_sch - B_sch/D_sch*C_sch);
            
            tmp2 = Trpr(free_dofs);
            
            [x, is_conv, rel_res, n_iters] = ...
                pcg(tmp, z, 1e-10, 100, model.L, model.U, tmp2(fem_dofs));
            
            y = D_sch \ (Fb_sch - C_sch * x);
            
            Tr(free_dofs) = [x; y];
        else
            error('Unknown solver type')
        end
        
        res = norm(Tr - Trpr);
        if (res < max_res)
            break
        end
        
        
    end
else
    
end
end
