function [P11, P12, P21, P22] = assembleContactFemFvm(local_contact, T_prev)
coder.allowpcode('plain')

th_body_dofs = local_contact.th_dofs;
th_body_massflow_dofs = local_contact.th_massflow_dofs;

temp_data = local_contact.th_data_T_ene;
rho_data = local_contact.th_data_T_rho;
nu_data = local_contact.th_data_T_nu;
mi_data = local_contact.th_data_T_mi;
lambda_data = local_contact.th_data_T_lambda;
cp_data = local_contact.th_data_T_cp;
th_A = local_contact.th_A;
th_radius = local_contact.th_radius;

fem_body_dofs = local_contact.fem_dofs;

energy = T_prev(th_body_dofs);

temp = interp1(temp_data(:,2), temp_data(:,1), energy);
cp = interp1(cp_data(:,1), cp_data(:,2), temp);

body_temp_full = local_contact.matrix_mapping_ur_to_pipe*T_prev(fem_body_dofs);

pipe_temp_full_projected = temp;

rho_pipe = interp1(rho_data(:,1), rho_data(:,2), temp);

velocity = T_prev(th_body_massflow_dofs) ... % TODO - fix, change
    ./ ([rho_pipe(1); (rho_pipe(1:end-1) + rho_pipe(2:end))./2; rho_pipe(end)]) ...
    ./ th_A;
mean_velocity = (velocity(1:end-1) + velocity(2:end))./2;
mean_velocity_projected = mean_velocity;

D_H = (th_radius(1:end-1) + th_radius(2:end)); % (r1+r2)/2 * 2
D_H_projected = D_H;

nu = interp1(nu_data(:,1), nu_data(:,2), pipe_temp_full_projected);

reynolds_number = abs(mean_velocity_projected) .* D_H_projected ./ nu;

% nusselt_number = self.computeNusseltNumber(reynolds_number, T_pipe_u_full_projected, body_u_full);
mi_wall = interp1(mi_data(:,1), mi_data(:,2), body_temp_full);
mi_fluid = interp1(mi_data(:,1), mi_data(:,2), pipe_temp_full_projected);

Pr_wall = interp1(cp_data(:,1), cp_data(:,2), body_temp_full) .* ...
    interp1(rho_data(:,1), rho_data(:,2), body_temp_full) .* ...
    interp1(nu_data(:,1), nu_data(:,2), pipe_temp_full_projected) ./ ...
    interp1(lambda_data(:,1), lambda_data(:,2), body_temp_full);

Pr_fluid = interp1(cp_data(:,1), cp_data(:,2), pipe_temp_full_projected) .* ...
    interp1(rho_data(:,1), rho_data(:,2), pipe_temp_full_projected) .* ...
    interp1(nu_data(:,1), nu_data(:,2), pipe_temp_full_projected) ./ ...
    interp1(lambda_data(:,1), lambda_data(:,2), pipe_temp_full_projected);

nusselt_number = zeros(size(reynolds_number));

laminar_ind = reynolds_number < 2e3;
transition_ind = reynolds_number >= 2e3 & reynolds_number < 1e4;
turbulent_ind =  reynolds_number > 1e4;

nusselt_number(laminar_ind) = 0.15 .* reynolds_number(laminar_ind).^0.33 .* ...
    Pr_fluid(laminar_ind).^0.43 .* (Pr_fluid(laminar_ind)./Pr_wall(laminar_ind)).^0.25;

Re_values = [1.9 2.2, 2.3, 2.5, 3.0, 3.5, 4, 5, 6, 7, 8, 9, 10] .* 1000;
K0_value = [1.9 2.2, 3.6, 4.9, 7.5, 10, 12.2, 16.5, 20, 24, 27, 30, 33];

nusselt_number(transition_ind) = interp1(Re_values,K0_value,reynolds_number(transition_ind)) .* ...
    Pr_fluid(transition_ind).^0.43 .* (Pr_fluid(transition_ind)./Pr_wall(transition_ind)).^0.25;

nusselt_number(turbulent_ind) = 0.023 .* reynolds_number(turbulent_ind).^0.8 .* ...
    Pr_fluid(turbulent_ind).^0.33 .* (mi_fluid(turbulent_ind)./mi_wall(turbulent_ind)).^0.14;

alpha = nusselt_number .* interp1(lambda_data(:,1), lambda_data(:,2), pipe_temp_full_projected) ./ D_H_projected;

alpha_projected = local_contact.projection_operator_pipe * alpha;
number_of_modes = size(local_contact.reduced_matrix_P11_basis, 1);

number_of_fv = size(th_A, 1) - 1; % number of volumes

% Init
reduced_matrix_P11 = zeros(number_of_modes, number_of_modes);
reduced_matrix_P12 = zeros(number_of_modes, number_of_fv);
reduced_matrix_P21 = zeros(number_of_fv, number_of_modes);
reduced_matrix_P22 = zeros(number_of_fv, number_of_fv);

coder.unroll(false)
for o = 1:numel(alpha_projected)
    reduced_matrix_P11 = reduced_matrix_P11 + local_contact.reduced_matrix_P11_basis(:, :, o)...
        *alpha_projected(o);
    reduced_matrix_P21 = reduced_matrix_P21 + local_contact.reduced_matrix_P21_basis(:, :, o)...
        *alpha_projected(o);
    reduced_matrix_P12 = reduced_matrix_P12 + local_contact.reduced_matrix_P12_basis(:, :, o)...
        *alpha_projected(o);
    reduced_matrix_P22 = reduced_matrix_P22 + local_contact.reduced_matrix_P22_basis(:, :, o)...
        *alpha_projected(o);
end

P11 = reduced_matrix_P11;
P12 = reduced_matrix_P12./repmat(cp', size(reduced_matrix_P12,1),1);
P21 = reduced_matrix_P21;
P22 = reduced_matrix_P22./repmat(cp', size(reduced_matrix_P22,1),1);

end