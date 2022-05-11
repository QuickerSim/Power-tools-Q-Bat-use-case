function [reduced_lambda, n_lambdas] = interpolateReducedLambda(fem_grid, modes_ids, reduced_u)
coder.allowpcode('plain')
sampled_u = fem_grid.sampling_matrix_diffusion_coeff_reduced(:, modes_ids) * reduced_u;
sampled_lambda = interploateLambda(fem_grid.lambda, sampled_u, fem_grid.flag_is_cylindrical);

if(isempty(fem_grid.recovery_matrix_diffusion_coeff))
    reduced_lambda = zeros(1, 1);
    n_lambdas = 1;
    return
end

if (size(sampled_lambda,1) == 3 &&  size(sampled_lambda,2) == 3)
    
    temp_lambda = [squeeze(sampled_lambda(1,1,:));...
        squeeze(sampled_lambda(2,2,:));...
        squeeze(sampled_lambda(3,3,:))];
    reduced_lambda = fem_grid.recovery_matrix_diffusion_coeff * temp_lambda;
    n_lambdas = 3;
    
elseif numel(sampled_lambda) > numel(sampled_u)
    
    reduced_lambda = fem_grid.recovery_matrix_diffusion_coeff * sampled_lambda(:);
    n_lambdas = 3;
    
else
    
    reduced_lambda = fem_grid.recovery_matrix_diffusion_coeff * sampled_lambda(:);
    n_lambdas = 1;
    
end

end