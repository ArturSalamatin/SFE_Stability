clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = -1;
sigma_min_limit = -3;
q = 1.003;
xBarLeft = 3e-3;
xBarRight = 1e-3;
%% packed bed params
% a0 = 0.2;
% alpha = [0, 0.1, 0.3, 0.5, 0.7];
% marker = ['o','s','d','*','x'];
%
% params = poly_case(a0, 0.3, (a0*a0/2)*1.1005);
% params.marker = marker(2);

alpha = 0.0;
a0 = 0.1;
tau0 = 0.4;
R = 2;
h = 0.001;
sigma = -1.821098878268054;

alpha = 0.5;
a0 = 0.3;
tau0 = 0.4;
R = 0.0001;
h = 5;
sigma = -1.62442043;

% alpha = 0.2;
% a0 = 0.01;
% tau0 = 0.47;
% R = 0.000001;
% h = 5;
% sigma = -2;

for Sigma = [sigma]
    sigma = Sigma;
    

params = poly_case(a0, alpha, tau0, R, h);
mesh = set_left_mesh(700, params, xBarLeft);

    (2+sigma)/params.C2
% solver = @(problem, mesh) solver_KellerBox(problem, mesh);
% solver = @(problem, mesh) solver_BVP(problem, mesh, params);

solver = @(problem, mesh) solver_RK(problem, mesh, params);
starter = @(sigma, params) starter_full(...
    solver, sigma, params, mesh);
sol = starter(sigma, params);
% sol.y(end,2)
pen = set_pen('r', '--');
plot_solution(pen, params, sol)
accuracy = (sol.condition/sol.y(sol.id,1)-1)

h_phi_div_gamma = sol.y(end, 3)*h/sol.y(end,4)

xi0 = params.z0/params.z2;
v = -params.g0/params.C1*(params.C2 + (1+sigma)*(1-xi0));
figure(701)
plot(xi0, v, 'd')
end
