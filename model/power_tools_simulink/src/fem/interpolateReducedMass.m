function reduced_mass_coeff = interpolateReducedMass(fem_grid, modes_ids, reduced_u)
coder.allowpcode('plain')
sampled_u = fem_grid.sampling_matrix_mass_coeff_reduced(:, modes_ids) * reduced_u;

if isscalar(fem_grid.rho.rho)
    sampled_rho = fem_grid.rho.rho * ones(size(sampled_u));
elseif length(fem_grid.rho.rho) > 1
    sampled_rho = nearestInterp1(fem_grid.rho.T, fem_grid.rho.rho, sampled_u);
else
    sampled_rho = NaN;
    error('non handled case');
end

if isscalar(fem_grid.cp.cp)
    sampled_cp = fem_grid.cp.cp * ones(size(sampled_u));
elseif length(fem_grid.cp.cp) > 1
    sampled_cp = nearestInterp1(fem_grid.cp.T, fem_grid.cp.cp, sampled_u);
else
    sampled_cp = NaN;
    error('non handled case');
end

sampled_mass = sampled_rho .* sampled_cp;
reduced_mass_coeff = fem_grid.recovery_matrix_mass_coeff * sampled_mass(:);

end