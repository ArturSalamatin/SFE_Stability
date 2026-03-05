% calculates eigenfunctions for a set of
% process parameters for the full problem, R > 0,
% of 4 searched-for functions
% using one of the solvers,
% BVP(from MatLab), RK(from MatLab), KellerBox(manually implemented)

clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = -1;
sigma_min_limit = -3;
q = 1.005;
xBarLeft = 1e-2;
xBarRight = 1e-4;
%% packed bed params
% a0 = 0.2;
% alpha = [0, 0.1, 0.3, 0.5, 0.7];
% marker = ['o','s','d','*','x'];
%
% params = poly_case(a0, 0.3, (a0*a0/2)*1.1005);
% params.marker = marker(2);

% alpha = 0.0;
% a0 = 0.1;
% tau0 = 0.4;
% R = 2;
% h = 0.001;
% sigma = -1.821098878268054;

alpha = 0.5;
a0 = 0.3;
tau0 = 0.4;
R = 0.0000000001;
h = 5;
sigma = -1.62442043;

% alpha = 0.1;
% a0 = 0.1;
% tau0 = 0.4;
% R = 0;
% sigma = -2.055860009887601;

% alpha = 0.2;
% a0 = 0.01;
% tau0 = 0.47;
% R = 0.000001;
% h = 5;
% sigma = -2;
N = 4000;
for Sigma = [sigma]
    %% set parameters
    sigma = Sigma;
    params = poly_case(a0, alpha, tau0, R, h);
    xi0 = params.z0/params.z2;
    v = -params.g0/params.C1*(params.C2 + (1+sigma)*(1-xi0));
    figure(701)
    hold on
    plot(xi0, v, 'd')
    x_left = linspace(params.a, params.a0, 10000);
    %% solve by KellerBox method
    disp(' ');
    disp('Keller Box solver:');
    mesh = set_full_mesh(N, params, 0);
    solver = @(problem, mesh) solver_KellerBox(problem, mesh, params);
    % solve the problem
    starter = @(sigma, params) starter_full(...
        solver, sigma, params, mesh);
    sol = starter(sigma, params);
    % plot the problem solution
    pen = set_pen('r', '--');
    plot_solution(pen, params, sol, x_left, 4)
    accuracy = (sol.condition/sol.y(sol.id,1)-1);
    disp(['accuracy = ', num2str(accuracy)]);
    %% solve by RK method
    disp('RK solver:');
    mesh = set_full_mesh(N, params, xBarLeft);
    solver = @(problem, mesh) solver_RK(problem, mesh, params);
    % solve the problem
    starter = @(sigma, params) starter_full(...
        solver, sigma, params, mesh);
    sol = starter(sigma, params);
    % plot the problem solution
    pen = set_pen('k', '--');
    plot_solution(pen, params, sol)
    accuracy = (sol.condition/sol.y(sol.id,1)-1);
    disp(['accuracy = ', num2str(accuracy)]);
    
    %% solve by BVP method
%     disp(' ');
%     disp('RK solver:');
%     mesh = set_full_mesh(N, params, xBarLeft);
%     solver = @(problem, mesh) solver_BVP(problem, mesh, params);
%     % solve the problem
%     starter = @(sigma, params) starter_full(...
%         solver, sigma, params, mesh);
%     sol = starter(sigma, params);
%     % plot the problem solution
%     pen = set_pen('k', '--');
%     plot_solution(pen, params, sol)
%     accuracy = (sol.condition/sol.y(sol.id,1)-1);
%     disp(['accuracy = ', num2str(accuracy)]);
end
