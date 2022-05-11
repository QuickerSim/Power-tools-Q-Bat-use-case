classdef simulink_fem_model < matlab.System & matlab.system.mixin.Propagates
    
    properties(Constant)
        % 10-Feb-2022 10:14:48
        
        % Number of solid components
        fem_n_ele = 28;
        
        % Number of cell components
        cell_n_ele = 5;
        
        % Number of pipe components
        th_n_ele = 0;
        
        % Number of air components
        fv_n_ele = 0;
        % Number of global dofs
        number_global_dofs = 1693;
    end
    
    properties(Nontunable)
        T_init_fem = NaN;
        T_init_th = NaN;
        T_init_fv = NaN;
        dt = NaN;
        max_res = NaN;
        max_sub_iter = NaN;
        n_fv_gradient_subiters = NaN;
        solver_type = NaN;
    end
    
    % Pre-computed constants
    properties(Access = private)
        model
        T_prev
        total_time = 0;
        fv_gradient;
    end
    
    properties(Constant, Access = private)
        fem_max_temp_to_save = [12, 13, 14, 15, 28];
        fem_mean_temp_to_save = [12, 13, 14, 15, 28];
        fem_min_temp_to_save = [12, 13, 14, 15, 28];
        fem_probed_temp_to_save = [];
        th_outlet_temp_to_save = [];
        fv_max_temp_to_save = [];
        fv_mean_temp_to_save = [];
        fv_min_temp_to_save = [];
        fv_probed_temp_to_save = [];
        fv_outlet_temp_to_save = [];
    end
    
    methods(Access = public)
        function obj = simulink_fem_model()
        end
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            data = load("data/fem_model.mat");
            
            obj.model = data.fem_model;
        end
        
        function [fem_mean_temp , fem_max_temp_output , fem_mean_temp_output , fem_min_temp_output       ] = ...
                stepImpl(obj, ele_source_term, bus_input)
            
            [Tr, n_fem_bodies, n_fem_grids, n_th_bodies, n_fv_grids, n_fv_bodies, fem_body_names, ...
                fem_grid_names, th_body_names, fv_grid_names, fv_body_names, fv_gradient] = doStep(obj.solver_type, ...
                obj.dt, obj.max_sub_iter, obj.T_prev, obj.total_time, obj.max_res, obj.model, obj.number_global_dofs, ...
                obj.n_fv_gradient_subiters, obj.fv_gradient, ele_source_term, bus_input);
            
            obj.T_prev = Tr;
            obj.fv_gradient = fv_gradient;
            
            fem_mean_temp = zeros(obj.cell_n_ele, 1);
            fem_max_temp_output = zeros(5,1);
            fem_min_temp_output = zeros(5,1);
            fem_mean_temp_output = zeros(5,1);
            fem_probed_temp_output = zeros(0,1);
            th_outlet_temp_output = zeros(0,1);
            fv_max_temp_output = zeros(0,1);
            fv_min_temp_output = zeros(0,1);
            fv_mean_temp_output = zeros(0,1);
            fv_probed_temp_output = zeros(0,1);
            fv_outlet_temp_output = zeros(0,1);
            
            
            
            
            
            for i = 1 : n_fem_bodies
                g_id = obj.model.fem_bodies.(fem_body_names{i}).grid_name;
                for j = 1 : n_fem_grids
                    if(j == g_id)
                    if any(i == unique([obj.fem_max_temp_to_save, obj.fem_min_temp_to_save]) ) 
    modes_ids = obj.model.fem_bodies.(fem_body_names{i}).modes_ids; 
    u_full = obj.model.fem_grids.(fem_grid_names{j}).reduced_basis(:, modes_ids) * ... 
        obj.T_prev(obj.model.fem_bodies.(fem_body_names{i}).lrgr);
    fem_max_temp_output(i == obj.fem_max_temp_to_save) = max(u_full);
    fem_min_temp_output(i == obj.fem_min_temp_to_save) = min(u_full);
end

                    
                    if ~isempty(obj.model.fem_bodies.(fem_body_names{i}).ele_body_id) || any(i == obj.fem_mean_temp_to_save)
                       average_operator =...
                          obj.model.fem_grids.(fem_grid_names{j}).average_operator(obj.model.fem_bodies.(fem_body_names{i}).modes_ids);
                       if ~isempty(obj.model.fem_bodies.(fem_body_names{i}).ele_body_id)
                          fem_mean_temp(obj.model.fem_bodies.(fem_body_names{i}).ele_body_id) = average_operator * obj.T_prev(obj.model.fem_bodies.(fem_body_names{i}).lrgr);
                       end
                       if any(i == obj.fem_mean_temp_to_save)
    fem_mean_temp_output(i == obj.fem_mean_temp_to_save) = average_operator * obj.T_prev(obj.model.fem_bodies.(fem_body_names{i}).lrgr);
end

                       average_operator = [];
                    end
                
                    end
                end
            end
            
            
            
            
            
            obj.total_time = obj.total_time + obj.dt;
        end
        
        function resetImpl(obj)
            obj.T_prev = doReset(obj.number_global_dofs, obj.model, obj.T_init_fem, obj.T_init_th, obj.T_init_fv);
            if obj.fv_n_ele > 0
                names = fieldnames(obj.model.fv_bodies);
                obj.fv_gradient = cell(1, obj.fv_n_ele);
                for i = 1 : obj.fv_n_ele
                    obj.fv_gradient{i} = zeros(3, numel(obj.model.fv_bodies.(names{i}).lrgr));
                end
            else
                obj.fv_gradient = 0;
            end
        end
        
        function num = getNumOutputsImpl(~)
            num = 4;
        end
        
        function [dataout , dataout_max , dataout_mean , dataout_min       ] = getOutputDataTypeImpl(~)
            dataout = 'double';
            dataout_max = 'double';
            dataout_mean = 'double';
            dataout_min = 'double';
            
            
            
            
            
            
            
        end
        
        function [cplxout , cplx_max , cplx_mean , cplx_min       ] = isOutputComplexImpl(~)
            cplxout = false;
            cplx_max = false;
            cplx_mean = false;
            cplx_min = false;
            
            
            
            
            
            
            
        end
        
        function [sz1 , sz_max , sz_mean , sz_min       ] = getOutputSizeImpl(obj)
            sz1 = [obj.cell_n_ele, 1];
            sz_max = [5, 1];
            sz_mean = [5, 1];
            sz_min = [5, 1];
            
            
            
            
            
            
            
        end
        
        function [fz1 , fz_max , fz_mean , fz_min       ] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz_max = true;
            fz_mean = true;
            fz_min = true;
            
            
            
            
            
            
            
        end
        
    end
    
    methods(Access = private)
    end
end