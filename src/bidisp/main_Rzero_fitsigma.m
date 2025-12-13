clc
clear all
% close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = 3;
sigma_min_limit = -2;
q = 1.008;
xBarLeft = 5e-3;
xBarRight = 0*4e-2;


sigma_guess = 1;
alpha = 0.5;
a0 = 0.2;
tau0 = 0.3;
R = 0;

sigma = 1.0889;
alpha = 0.1;
a0 = 0.1;
tau0 = 0.1;
R = 0;

alpha = 0.5;
a0 = 0.2;
tau0 = 0.35;
R = 0;
sigma = -0.124996446879085; % if plus

alpha = 0.18;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.264900506177209;

alpha = 0.17;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.317466186083304;

alpha = 0.1;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.717977418539169;


alpha = 0.2;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.162932030910141;

alpha = 0.2;
a0 = 0.1;
tau0 = 0.47;
R = 0;
sigma = -2.031115;
sigma = 2.151615819264947;

alpha = 0.2;
a0 = 0.1;
tau0 = 0.2864;
R = 0;
sigma_guess = 1.458503275213755;

alpha = 0.2;
a0 = 0.05;
tau0 = 0.2864;
R = 0;
sigma_guess = 2.276089956284676;

alpha = 0.2;
a0 = 0.03;
tau0 = 0.2864;
R = 0;
sigma_guess = 2.826374592543163;

alpha = 0.2;
a0 = 0.01;
tau0 = 0.2864;
R = 0;
sigma_guess = 3.652738419309777;

params = poly_case(a0, alpha, tau0, R);
params.pen = set_pen('m', '-');
mesh = set_left_mesh(100, params, xBarLeft);
solver = @(problem, mesh) solver_RK(problem, mesh, params);
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
% for i = 2:numel(sigma)
%     if(out(i) < out(i-1))
%         out(i-1) = NaN;
%         break;
%     end
% end
%% plot F(sigma)
figure(3000)
hold on
axis([sigma_min_limit sigma_max_limit -1 1])
plot(sigma, out, 'b-', 'LineWidth', 1)
hold on
grid on
%% fit sigma
[sigma, sol] = fit_sigma(starter, params, sigma_guess);

plot_solution(params.pen, params, sol)