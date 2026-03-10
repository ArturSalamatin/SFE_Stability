clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = 3;
sigma_min_limit = -2.1;
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
sigma_guess = -2.031672191741087;

alpha = 0.2;
a0 = 0.01;
tau0 = 0.47;
R = 0;
sigma_guess = -2.031672191741087;
% sigma_guess = -15.588081024788260;

% alpha = 0.2;
% a0 = 0.1;
% tau0 = 0.2864;
% R = 0;
% sigma_guess = -1.990654452314103;

% alpha = 0.2;
% a0 = 0.05;
% tau0 = 0.2864;
% R = 0;
% sigma_guess = -2.017528213766199;
%%
alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 0;
h = 5;
sigma_guess = -2.055860009887601;


alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 1;
h = 5;
sigma_guess = -1.301938994100083;

alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 1;
h = 50;
sigma_guess = -0.487594835246630;

alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 1;
h = 500;
sigma_guess = -0.151878663989612;

alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 1;
h = 1500;
sigma_guess = -0.125215532909594;

% alpha = 0.1;
% a0 = 0.1;
% tau0 = 0.4;
% R = 2;
% h = 5;
% sigma_guess = -0.245363645840584;

% alpha = 0.1;
% a0 = 0.1;
% tau0 = 0.4;
% R = 2;
% h = 65;
% sigma_guess = 1.794832507640907;


alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 1;
h = 500;
sigma_guess = -0.151878663989612;

params = poly_case(a0, alpha, tau0, R, h);
params.pen = set_pen('k', '-');
mesh = set_full_mesh(3500, params, 0);
solver = @(problem, mesh) solver_KellerBox(problem, mesh, params);
% solver = @(problem, mesh) solver_RK(problem, mesh, params);
% solver = @(problem, mesh) solver_BVP(problem, mesh, params);
starter = @(sigma, params) starter_full(...
    solver, sigma, params, mesh);
%% plot functional
% sigma = linspace(-3,1,81);
% out = zeros(size(sigma));
% for i = 1:numel(sigma)
%     sol = starter(sigma(i), params);
%     out(i) = sol.condition;
% end
% %% do not plot jumps
% for i = 2:numel(sigma)
%     if(out(i) < out(i-1))
%         out(i-1) = NaN;
%         break;
%     end
% end
% %% plot F(sigma)
% figure(3000)
% hold on
% % axis([sigma_min_limit sigma_max_limit -1 1])
% plot(sigma, out, 'r-', 'LineWidth', 1)
% hold on
% grid on
%% fit sigma
[sigma, sol] = fit_sigma(starter, params, sigma_guess);

x_left = linspace(params.a, params.a0, 10001);
plot_solution(params.pen, params, sol, x_left, 4)
accuracy = (sol.condition/sol.y(sol.id,1)-1);
disp(['accuracy = ', num2str(accuracy)]);

if(params.R < 1e-7)
    xi0 = params.z0/params.z2;
    v = -params.g0/params.C1*(params.C2 + (1+sigma)*(1-xi0));
    figure(701)
    plot(xi0, v, 'd')
end