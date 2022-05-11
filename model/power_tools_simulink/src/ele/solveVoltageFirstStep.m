function u_voltage = solveVoltageFirstStep(ele_bodies, u_ele_old, temperature)
coder.allowpcode('plain')
names = fieldnames(ele_bodies);
n_electro_bodies = numel(names);
u_voltage = zeros(numel(u_ele_old)/2, 1);

for i = 1 : n_electro_bodies
    curr_ele_body = ele_bodies.(names{i});
    
    lrgr = curr_ele_body.lrgr;
    curr_temp = temperature(curr_ele_body.fem_body_id);
    SOC = u_ele_old(lrgr(2));
    
    if numel(curr_ele_body.U_OCV.U_OCV ) == 1
        u_voltage(i) = curr_ele_body.U_OCV.U_OCV;
    else
        u_voltage(i) = nearestInterp2(curr_ele_body.U_OCV.T, ...
            curr_ele_body.U_OCV.SOC, ...
            curr_ele_body.U_OCV.U_OCV, ...
            curr_temp, ...
            SOC);
    end
end

end