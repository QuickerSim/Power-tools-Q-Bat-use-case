classdef simulink_electro_model < matlab.System & matlab.system.mixin.Propagates
    
    properties(Constant)
        % 10-Feb-2022 10:14:48
        
        % Number of electro components
        n_ele = 5;
        
        active = 1
    end
    
    properties(Nontunable)
        dt = NaN;
        n_substeps = NaN;
        SOC_init = NaN;
    end
    
    properties(DiscreteState)
        u_ele
        u_ele_old
    end
    
    % Pre-computed constants
    properties(Access = private)
        
        model
        source
        source_term_aggregated
        prev_source_term
        curr_substep
        total_time = 0
        
    end
    
    properties(Constant, Access = private)
        electro_SOC_to_save = [1, 2, 3, 4, 5];
        electro_vol_to_save = [1];
    end
    
    methods(Access = public)
        function obj = simulink_electro_model()
        end
    end
    
    methods(Access = protected)
        
        function setupImpl(obj)
            data = load("data/electro_model.mat");
            
            obj.model = data.electro_model;
        end
        
        function [source_term , electro_SOC_output , electro_vol_output] = stepImpl(obj, current, temperature)
            if obj.active == 0 % case when there is no electro bodies
                source_term = zeros(obj.n_ele, 1);
                obj.u_ele = obj.u_ele + 1;
                obj.u_ele_old = obj.u_ele_old + 1;
            else
                
                first_simulink_cycle = false;
                if obj.total_time == 0
                    first_simulink_cycle = true;
                end
                
                if ~first_simulink_cycle
                    obj.source = zeros(obj.n_ele, 1);
                    
                    [u_electro, u_voltage] = obj.solveElectroSolution(obj.dt, current, temperature);
                    obj.source_term_aggregated = obj.source_term_aggregated + ...
                        obj.source;
                    
                    if(obj.curr_substep == obj.n_substeps)
                        mean_source = obj.source_term_aggregated./obj.n_substeps;
                        obj.prev_source_term = mean_source;
                        obj.source_term_aggregated = zeros(obj.n_ele, 1);
                        obj.curr_substep = 1;
                    else
                        mean_source = zeros(obj.n_ele, 1);
                        obj.curr_substep = obj.curr_substep + 1;
                    end
                    
                else 
                    mean_source = zeros(obj.n_ele, 1);
                    u_voltage = obj.solveVoltageFirstStep(temperature);
                end
                
                source_term = mean_source;
                
                electro_names = fieldnames(obj.model.ele_bodies);
                n_electro = numel(electro_names);
                
                electro_SOC_output = zeros(5,1);
                electro_vol_output = zeros(1,1);
                
                for i = 1 : numel(electro_names)
if any(i == obj.electro_SOC_to_save)
    local_electro = obj.model.ele_bodies.(electro_names{i});
     local_dofs = local_electro.lrgr;
      electro_SOC_output(i == obj.electro_SOC_to_save) = obj.u_ele(local_dofs(2));
   end
end
                for branch_id = 1 : numel(electro_vol_output)
     electro_vol_output(branch_id) = sum(u_voltage(obj.model.ele_contacts.electro_contacts_1.vol_lrgr(1, :)));
end
                
                    
                obj.total_time = obj.total_time + obj.dt;
            end
            
        end
        
        function resetImpl(obj)
            
            obj.u_ele = zeros(2*obj.n_ele, 1);
            obj.u_ele_old = zeros(2*obj.n_ele, 1);
            obj.u_ele(2:2:end) = obj.SOC_init;
            obj.u_ele_old(2:2:end) = obj.SOC_init;
            obj.source = zeros(obj.n_ele, 1);
            
            obj.source_term_aggregated = zeros(obj.n_ele, 1);
            obj.prev_source_term = zeros(obj.n_ele, 1);
            obj.curr_substep = 1;
        end
        
        function num = getNumOutputsImpl(~)
            num = 3;
        end
        
        function [dataout , dataout_SOC , dataout_vol] = getOutputDataTypeImpl(~)
            dataout = 'double';
            dataout_SOC = 'double';
            dataout_vol = 'double';
        end
        
        function [cplxout , cplx_SOC , cplx_vol] = isOutputComplexImpl(~)
            cplxout = false;
            cplx_SOC = false;
            cplx_vol = false;
        end
        
        function [sz,dt,cp] = getDiscreteStateSpecificationImpl(obj, name)
            switch name
                case {'u_ele', 'u_ele_old'}
                    sz = [2*obj.n_ele, 1];
                    dt = "double";
                    cp = false;
                    
                case 'source'
                    sz = [obj.n_ele, 1];
                    dt = "double";
                    cp = false;
                    
                otherwise
                    sz = [1, 1];
                    dt = "double";
                    cp = false;
            end
            
        end
        
        function [sz1 , sz_SOC , sz_vol] = getOutputSizeImpl(obj)
            sz1 = [obj.n_ele];
            sz_SOC = [5, 1];
            sz_vol = [1, 1];
        end
        
        function [fz1 , fz_SOC , fz_vol] = isOutputFixedSizeImpl(~)
            fz1 = true;
            fz_SOC = true;
            fz_vol = true;
        end
        
    end
    
    
    methods(Access = private)
        
        function [u_electro, u_voltage] = solveElectroSolution(self, dt, current, temperature)
            [u_electro, u_voltage, self.source, self.u_ele, self.u_ele_old] = ...
                solveElectroSolution(self.model.ele_bodies, self.u_ele_old, self.u_ele, self.source, dt, current, temperature);
        end
        
       function u_voltage = solveVoltageFirstStep(self, temperature)
             u_voltage = solveVoltageFirstStep(self.model.ele_bodies, self.u_ele_old, temperature);
        end
        
    end
    
end