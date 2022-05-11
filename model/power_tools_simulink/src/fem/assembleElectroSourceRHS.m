function F_ele = assembleElectroSourceRHS(fem_grids, fem_bodies, number_global_dofs, source_values)
coder.allowpcode('plain')
body_names = fieldnames(fem_bodies);
n_fem_bodies = numel(body_names);

grid_names = fieldnames(fem_grids);
n_fem_grids = numel(grid_names);
F_ele = zeros(number_global_dofs, 1);

for i = 1 : n_fem_bodies
    current_body = fem_bodies.(body_names{i});
    
    if(~isempty(current_body.F_ele))
        
        g_id = current_body.grid_name;
        
        for j = 1 : n_fem_grids
            if(j == g_id)
                curr_volume = fem_grids.(grid_names{j}).grid_volume;
                curr_ele_body_id = current_body.ele_body_id;
                curr_lrgr = current_body.lrgr;
                curr_Fr = current_body.F_ele;
                source_values(curr_ele_body_id) = ...
                    source_values(curr_ele_body_id)/curr_volume;
                F_ele(curr_lrgr) = F_ele(curr_lrgr) + ...
                    curr_Fr*source_values(curr_ele_body_id);
            end
        end
        
    end
    
end

end