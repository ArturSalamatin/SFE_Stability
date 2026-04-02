clc
% close all
clear all

a0 = 0.2;
R = 1;
alpha = 0.5;
B = 1;

params = poly_case(a0, alpha, R, B);


params.L = 1.1; % spatial dimension
params.z_step = 0.0005;%0.001/2/2; % spatial grid step size
params.alpha = 1; % Bessel function eigenvalue
params.omega = 2; % velocity perturbation frequency
params.T = 0.5;
params.period = 0.5;
params.Nt = 41;
params.color = 'b-';



params.t_nodes_count = 1 + ...
    round(params.omega*(params.Nt-1)*params.T/params.period);
params.t_grid = linspace(1e-5,params.T, params.t_nodes_count);

params.t_plot_nodes = 1:(params.Nt-1)/4:params.t_nodes_count;

[params.z_step, params.z] = z_grid(params.z_step, params.L);
params.m = numel(params.z); % nmbr of z-nodes
[u0, du0] = init_value(params.z, params.omega);
params.dpsi0 = du0(params.m+1:end);
params.history = zeros(params.m, 1);

out = IVBP_solver(u0, params);

z = params.z;
figure(101)
hold on
axis([0 3 -Inf Inf])
plot(z, out(1:params.m,params.t_plot_nodes), params.color)

figure(201)
hold on
axis([0 1 -Inf Inf])
plot(z, out(params.m+1:end,params.t_plot_nodes), params.color)


t = params.t_grid;
m = params.m;
figure(301)
hold on
axis([0 0.5 -Inf Inf])
plot(t, max(abs(out(m+(1:m), :)), [], 1), params.color)

% psi = out(m+(1:m), :);
% psi_mid = (psi(:,1:end-1)+psi(:,2:end))/2;
% tau = diff(t);
% t_mid = (t(2:end)+t(1:end-1))/2;
% dpsi = diff(psi, 1, 2);
% dpsidt = dpsi./repmat(tau, params.m, 1);
% 
% 
% sigma_t = dpsidt./psi_mid;
% sigma_sqrt = sqrt(sum(sigma_t.^2, 1, 'omitnan')./sum(~isnan(sigma_t), 1));
% figure(401)
% hold on
% plot(t_mid, sigma_sqrt, params.color)





