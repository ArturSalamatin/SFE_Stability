clc
clear all
% close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = -1;
sigma_min_limit = -3;
q = 1.008;
xBarLeft = 5e-3;
xBarRight = 0*4e-2;


sigma_guess = -1.604021809231366;
alpha = 0.5;
a0 = 0.2;
tau0 = 0.3;
R = 0;

sigma = -2.05;
% sigma = -1.604021809231366;
% sigma = -1.9999;
alpha = 0.1;
a0 = 0.1;
tau0 = 0.1;
R = 0;

params = poly_case(a0, alpha, tau0, R);
params.pen = set_pen('b', '-');
mesh = set_left_mesh(2500, params, xBarLeft);
solver = @(problem) solver_KellerBox(problem, mesh);
starter = @(sigma, params) starter_R_zero_X_Psi(...
                solver, sigma, params, mesh);


%% plot functional
sigma = linspace(sigma_min_limit,sigma_max_limit,81);
out = zeros(size(sigma));
for i = 1:numel(sigma)
    sol = starter(sigma(i), params);
    out(i) = sol.condition/sol.y(end,1)-1;
end
%% do not plot jumps
for i = 2:numel(sigma)
    if(out(i) < out(i-1))
        out(i-1) = NaN;
        break;
    end
end
%% plot F(sigma)
figure(3000)
hold on
axis([sigma_min_limit sigma_max_limit -1 1])
plot(sigma, out, 'b-', 'LineWidth', 1)
hold on
grid on
%% fit sigma
[sigma, sol] = fit_sigma(starter, params, sigma_guess);

bvp4c
plot_solution(params.pen, params, sol)