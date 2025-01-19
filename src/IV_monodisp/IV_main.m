clc
% close all
clear all

params.L = 1.2; % spatial dimension
params.z_step = 0.001/2/2/2; % spatial grid step size
params.alpha = 100; % Bessel function eigenvalue
params.R = 0.001; % viscosity increment parameter
params.omega = 20; % perturbation frequency


% t = 0; % initial time moment

[params.z_step, params.z] = z_grid(params.z_step, params.L);
params.m = numel(params.z);
[u0, du0] = init_value(params.z, params.omega);
params.chi2 = u0(params.m+1:end);
params.dpsi0 = 0*du0(params.m+1:end);
params.history = zeros(params.m, 1);

out = IVBP_solver(0*u0, params);

diff(out, 2);



figure(101)
hold on
plot(params.z, out(1:params.m,1:200:end), 'b--')

figure(201)
hold on
plot(params.z, out(params.m+1:end,1:200:end), 'b--')


