function [M, C] = assembleTHBody(th_body, T_prev)
coder.allowpcode('plain')

lrgr_T = th_body.lrgr_T;
lrgr_massflow = th_body.lrgr_massflow;

energy = T_prev(lrgr_T);
mass_flow = T_prev(lrgr_massflow);

temp_data = th_body.th_pipe.data_T_ene;
temp = interp1(temp_data(:,2), temp_data(:,1), energy);
rho_data = th_body.th_pipe.data_T_rho;
rho = interp1(rho_data(:,1), rho_data(:,2), temp);

A = th_body.th_pipe.A;
V = th_body.th_pipe.V;
n_volumes = numel(A)-1;

% calculate velocities
vel = mass_flow ...
    ./ ([rho(1); (rho(1:end-1) + rho(2:end))./2; rho(end)]) ...
    ./ A;

% mass matrix
M = zeros(n_volumes, n_volumes);

% convection matrix
C = zeros(n_volumes, n_volumes);

M(1, 1) = rho(1) * V(1);
C(1, 1) = (rho(1)+rho(2))./2 * A(2) * vel(2);

coder.unroll(false)
for i = 2:n_volumes-1
    M(i, i) = rho(i) * V(i);
    C(i, i) = (rho(i)+rho(i+1)) ./2 * A(i+1) * vel(i+1);

    % M(i, i-1) = 0;
    C(i, i-1) = -(rho(i-1)+rho(i)) ./2 * A(i) * vel(i);
end

i = n_volumes;
M(i, i) = rho(i) * V(i);
C(i, i) = rho(i) * A(i+1) * vel(i+1);

% M(i, i-1) = 0;
C(i, i-1) = -(rho(i-1)+rho(i)) ./2 * A(i) * vel(i);


end