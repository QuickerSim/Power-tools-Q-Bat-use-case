function [M, K] = getNonstaticMatrices(grid, body, T_body_reduced)
coder.allowpcode('plain')
n = numel(body.lrgr);

M = zeros(n, n);
K = zeros(n, n);

if grid.flag_diffusion_non_linear

    [lambda_reduced, n_lambdas] = interpolateReducedLambda(grid, body.modes_ids, T_body_reduced);

    nKShapes = size(grid.reduced_matrix_basis_diffusion, 3);
    if (n_lambdas == 3)
        nKShapes = nKShapes / n_lambdas;
    end

    for lSId = 1 : nKShapes * n_lambdas
        K = K + ...
            grid.reduced_matrix_basis_diffusion(...
            body.modes_ids, body.modes_ids, lSId) ...
            * lambda_reduced(lSId);
    end

end

if grid.flag_mass_non_linear

    rho_cp_reduced = interpolateReducedMass(grid, body.modes_ids, T_body_reduced);
    num_of_something = size(grid.reduced_matrix_basis_mass, 3); % todo change this name

    if length(rho_cp_reduced) == num_of_something
        for lSId = 1 : num_of_something
            M = M + grid.reduced_matrix_basis_mass(body.modes_ids, body.modes_ids, lSId) ...
                * rho_cp_reduced(lSId);
        end
    end

end

end