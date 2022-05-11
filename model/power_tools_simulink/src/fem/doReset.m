function T_prev = doReset(number_global_dofs, model, T_init_fem, T_init_th, T_init_fv)
coder.allowpcode('plain')
% Assume uniform temperature distribution in whole battery
T_prev = zeros(number_global_dofs, 1);

fem_body_names = fieldnames(model.fem_bodies);
n_fem_bodies = numel(fem_body_names);

if isfield(model, "fv_bodies")
    fv_body_names = fieldnames(model.fv_bodies);
    n_fv_bodies = numel(fv_body_names);
else
    n_fv_bodies = 0;
end

if isfield(model, "th_bodies")
    th_body_names = fieldnames(model.th_bodies);
    n_th_bodies = numel(th_body_names);
else
    n_th_bodies = 0;
end

if(numel(T_init_fem) == 1)
    for i = 1 : n_fem_bodies
        local_dofs = model.fem_bodies.(fem_body_names{i}).lrgr;
        T_prev(local_dofs) = model.fem_bodies.(fem_body_names{i}).ur0 * T_init_fem;
    end
    
elseif(numel(T_init_fem) == n_fem_bodies)
    for i = 1 : n_fem_bodies
        local_dofs = model.fem_bodies.(fem_body_names{i}).lrgr;
        T_prev(local_dofs) = model.fem_bodies.(fem_body_names{i}).ur0 * T_init_fem(i);
    end
else
    error('Invalid T_init_fem size.')
end

if(numel(T_init_fv) == 1)
    for i = 1 : n_fv_bodies
        local_dofs = model.fv_bodies.(fv_body_names{i}).lrgr;
        T_prev(local_dofs) = T_init_fv;
    end
    
elseif(numel(T_init_fv) == n_fv_bodies)
    for i = 1 : n_fv_bodies
        local_dofs = model.fv_bodies.(fv_body_names{i}).lrgr;
        T_prev(local_dofs) = T_init_fv(i);
    end
elseif(n_fv_bodies > 0)
    error('Invalid T_init_fv size.')
end

for i = 1:n_th_bodies
    th_pipe = model.th_bodies.(th_body_names{i}).th_pipe;
    
    if numel(T_init_th) == 1
        temp_init = T_init_th;

        energy_data = th_pipe.data_T_ene;
        energy_init = interp1(energy_data(:,1), energy_data(:,2), temp_init);

        local_dofs_T = model.th_bodies.(th_body_names{i}).lrgr_T;
        local_dofs_massflow = model.th_bodies.(th_body_names{i}).lrgr_massflow;

        T_prev(local_dofs_T) = energy_init;
        T_prev(local_dofs_massflow) = model.th_bodies.(th_body_names{i}).mass_flow0;
    elseif numel(T_init_th) == n_th_bodies
        temp_init = T_init_th(i);

        energy_data = th_pipe.data_T_ene;
        energy_init = interp1(energy_data(:,1), energy_data(:,2), temp_init);

        local_dofs_T = model.th_bodies.(th_body_names{i}).lrgr_T;
        local_dofs_massflow = model.th_bodies.(th_body_names{i}).lrgr_massflow;

        T_prev(local_dofs_T) = energy_init;
        T_prev(local_dofs_massflow) = model.th_bodies.(th_body_names{i}).mass_flow0;
    else
        error('Invalid T_init_th size.')
    end
            
            
end
    
end