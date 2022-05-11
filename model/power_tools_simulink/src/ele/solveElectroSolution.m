function [u_electro, u_voltage, source, u_ele, u_ele_old] = solveElectroSolution(ele_bodies, u_ele_old, u_ele, source, dt, current, temperature)
coder.allowpcode('plain')

names = fieldnames(ele_bodies);
n_electro_bodies = numel(names);

u_voltage = zeros(numel(u_ele_old)/2, 1);

for i = 1 : n_electro_bodies
    curr_ele_body = ele_bodies.(names{i});
    model_type = curr_ele_body.model_type;
    lrgr = curr_ele_body.lrgr;
    curr_temp = temperature(curr_ele_body.fem_body_id);
    SOC = u_ele_old(lrgr(2));
    circuit_factor = curr_ele_body.circuit_factor;
    if numel(current) == 1
        curr_current = current * circuit_factor;
    else
        curr_current = current(i) * circuit_factor;
    end
    
    if numel(curr_ele_body.U_OCV.U_OCV ) == 1
        U_OCV = curr_ele_body.U_OCV.U_OCV;
    else
        U_OCV = nearestInterp2(curr_ele_body.U_OCV.T, ...
            curr_ele_body.U_OCV.SOC, ...
            curr_ele_body.U_OCV.U_OCV, ...
            curr_temp, ...
            SOC);
    end
    
    if(SOC <= 0) && (curr_current >= 0)
        u_ele(lrgr) = [0, 0];
        source(i) = 0;
        u_voltage(i) = 0;
        continue
        
    elseif(SOC > 1) && (curr_current <= 0)
        u_ele(lrgr) = [0, 1];
        source(i) = 0;
        u_voltage(i) = U_OCV;
        continue
        
    end
    
    A = curr_ele_body.A;
    
    if numel(curr_ele_body.R_1_charge.R_1) == 1 && numel(curr_ele_body.R_1_discharge.R_1) == 1
        R1 = (curr_current >= 0) * curr_ele_body.R_1_discharge.R_1 + ...
            (curr_current < 0) * curr_ele_body.R_1_charge.R_1;
    else
        R1 = (curr_current >= 0) * nearestInterp2(curr_ele_body.R_1_discharge.T, curr_ele_body.R_1_discharge.SOC, ...
            curr_ele_body.R_1_discharge.R_1, curr_temp, SOC) + ...
            (curr_current < 0) * nearestInterp2(curr_ele_body.R_1_charge.T, curr_ele_body.R_1_charge.SOC, ...
            curr_ele_body.R_1_charge.R_1, curr_temp, SOC);
    end
    
    
    switch model_type
        
        case 'R'
            R2 = 0;
            C = 0;
            
            M_loc = [0; -A];
            K_loc = [inf; 0];
            
            
        case 'RC'
            if numel(curr_ele_body.R_2_charge.R_2) == 1 && numel(curr_ele_body.R_2_discharge.R_2) == 1
                R2 = (curr_current >= 0) * curr_ele_body.R_2_discharge.R_2 + ...
                    (curr_current < 0) * curr_ele_body.R_2_charge.R_2;
            else
                R2 = (curr_current >= 0) * nearestInterp2(curr_ele_body.R_2_discharge.T, curr_ele_body.R_2_discharge.SOC, ...
                    curr_ele_body.R_2_discharge.R_2, curr_temp, SOC) + ...
                    (curr_current < 0) * nearestInterp2(curr_ele_body.R_2_charge.T, curr_ele_body.R_2_charge.SOC, ...
                    curr_ele_body.R_2_charge.R_2, curr_temp, SOC);
            end
            
            if numel(curr_ele_body.C_charge.C) == 1 && numel(curr_ele_body.C_discharge.C) == 1
                C = (curr_current >= 0) * curr_ele_body.C_discharge.C + ...
                    (curr_current < 0) * curr_ele_body.C_charge.C;
            else
                C = (curr_current >= 0) * nearestInterp2(curr_ele_body.C_discharge.T, curr_ele_body.C_discharge.SOC, ...
                    curr_ele_body.C_discharge.C, curr_temp, SOC) + ...
                    (curr_current < 0) * nearestInterp2(curr_ele_body.C_charge.T, curr_ele_body.C_charge.SOC, ...
                    curr_ele_body.C_charge.C, curr_temp, SOC);
            end
            
            M_loc = [1.0; -A];
            K_loc = [1.0/(R2*C); 0];
            
            
        otherwise
            error('unknown type of electric model');
    end
    
    u_ele(lrgr) = (M_loc.*u_ele_old(lrgr)/dt ...
        + curr_current ) ./ (M_loc/dt + K_loc);
    
    curr_ele_body.SOC = u_ele(lrgr(2));
    
    if (C == 0)
        source(i) = curr_current^2*R1;
        u_voltage(i) = U_OCV - curr_current*R1;
    else
        source(i) = curr_current^2*R1 + (u_ele(lrgr(1))/C)^2/R2;
        u_voltage(i) = U_OCV - curr_current*R1 - u_ele(lrgr(1))/C;
    end
    
end

u_ele_old = u_ele;
u_electro = u_ele;
end