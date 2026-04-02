clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = -1;
sigma_min_limit = -3;
q = 1.003;
xBarLeft = 3e-3;
xBarRight = 0;
%% packed bed params
% a0 = 0.2;
% alpha = [0, 0.1, 0.3, 0.5, 0.7];
% marker = ['o','s','d','*','x'];
%
% params = poly_case(a0, 0.3, (a0*a0/2)*1.1005);
% params.marker = marker(2);
%
% N = 150;
% % q = 1.06;
% % mesh = log_mesh(q, N, params);
% mesh = quasiuniform_mesh_Frobenius(1e-16, N, params);
% starter = @(sigma, params) starter_R_zero_X_Psi(...
%     sigma, params, mesh);

sigma = -2.02815;
sigma = -2.0553173;
sigma = -1.942; % a0 = 0.3
sigma = -1.98669; % a0 = 0.2; tau0 = 0.45;
sigma = -1.988;   % a0 = 0.2; tau0 = 0.47;
sigma = -2.097600; % a0 = 0.005; n = 100
sigma = -2.0920; % a0 = 0.005; n = 200
sigma = -2.0919; % a0 = 0.005; n = 2000
sigma = -1.986683;
%%
alpha = 0.2;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.162932030910141;
%%
%%
alpha = 0.1;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.717977418539169;
%%
%%
alpha = 0.17;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.317466186083304;
%%
%%
alpha = 0.18;
a0 = 0.2;
tau0 = 0.47;
R = 0;
sigma = 1.264900506177209;
%%
%%
% alpha = 0.2;
% a0 = 0.1;
% tau0 = 0.47;
% R = 0;
% sigma = -2.031115;
% sigma = 2.151615819264947;
%%
%%
% alpha = 0.5;
% a0 = 0.2;
% tau0 = 0.35;
% R = 0;
% sigma = -1.6338658; % if minus
% sigma = -0.124996446878894; % if plus
%%
alpha = 0.2;
a0 = 0.1;
tau0 = 0.2864;
R = 0;
sigma = -1.6338658;
sigma = -1.9906507;
%%

% close all

% alpha = 0.4;
% a0 = 0.1;
% tau0 = 0.4;
% R = 0;
% sigma = -1.821098878268054;
% % sigma = -15.588081024788260;

alpha = 0.2;
a0 = 0.01;
tau0 = 0.47;
R = 0;
sigma = -2.031672191741087;

alpha = 0.1;
a0 = 0.1;
tau0 = 0.4;
R = 0;
sigma = -2.055860009887601;

params = poly_case(a0, alpha, tau0, R);
mesh = set_left_mesh(700, params, xBarLeft);

% solver = @(problem, mesh) solver_KellerBox(problem, mesh);
solver = @(problem, mesh) solver_BVP(problem, mesh, params);
% starter = @(sigma, params) starter_R_zero_X_Psi(...
%     solver, sigma, params, mesh);
starter = @(sigma, params) starter_R_zero_X_Psi(...
    solver, sigma, params, mesh);
sol = starter(sigma, params);
% sol.y(end,2)
pen = set_pen('b', '--');
plot_solution(pen, params, sol)
accuracy = (sol.condition/sol.y(sol.id,1)-1)


solver = @(problem, mesh) solver_RK(problem, mesh, params);
starter = @(sigma, params) starter_R_zero_X_Psi(...
    solver, sigma, params, mesh);
sol = starter(sigma, params);
% sol.y(end,2)
pen = set_pen('k', '--');
plot_solution(pen, params, sol)


xi0 = params.z0/params.z2;
v = -params.g0/params.C1*(params.C2 + (1+sigma)*(1-xi0));
figure(701)
plot(xi0, v, 'd')

s=0;



% label = '';
% N = 801;
% for alpha = 0.2 % 0.2:0.2:0.8
%     num = 2200+alpha*10;
%     clc
%     col = {'k'};
%     style = {'-'};
%     pen = set_pen(col, style);
%     %% make calculation grid
%     eps = 10e-3;
%     a0 = 0.04;
%     tau0 = 0.35;
%     %     tau0 = linspace(eps, 0.5 - eps, 20);
%     %     a0 = linspace(eps,1-eps, 20);
%     [tau0, a0] = meshgrid(tau0, a0);
%     %% run calculations
%     [tau0, a0, sigma] = ...
%         calc_sigma_R_zero(N, alpha, tau0, a0, pen);
%     %     save(['R_zero_data/alpha_', num2str(alpha*100), '.mat'], ...
%     %         'tau0', 'a0', 'sigma', 'alpha')
%     %
%     %     figure(num)
%     %     hold on
%     %     axis([0 0.5 0 1])
%     %     contour(tau0, a0, sigma, ...
%     %         linspace(-2,-1,11),...
%     %         'ShowText','on',...
%     %         'linecolor', 'black')
%     %
%     %     x = linspace(0, 0.5, 700);
%     %     y = sqrt(2*x);
%     %     plot(x,y, 'k-', 'linewidth', 1)
%     %
%     %     x = -2:0.1:-1;
%     %     y = 1/(1-alpha)*(alpha./(x+2)-alpha);
%     %     z = y.*y/2;
%     %     plot(z,y,'ks', 'MarkerFaceColor', 'black')
%     %
%     %     switch alpha
%     %         case 0.2
%     %             label = 'A';
%     %         case 0.4
%     %             label = 'B';
%     %         case 0.6
%     %             label = 'C';
%     %         case 0.8
%     %             label = 'D';
%     %     end
%     %
%     %     box on
%     %     set(gca, 'xTick', [0:0.1:0.5])
%     %     set(gca, 'yTick', [0:0.2:1.0])
%     %     set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
%     %     set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])
%     %     ylabel(['{\ita}_0, ' char(8211)])
%     %     xlabel(['{\it\tau}_0, ' char(8211)])
%     %     annotation('textbox',...
%     %         [0.207441165123454 0.793200250350447 0.115740738211223 0.141129029253798],...
%     %         'String',{label},...
%     %         'LineStyle','none',...
%     %         'FontSize',16,...
%     %         'FontName','Times New Roman');
%     %
%     %     hFig = findobj('Type', 'figure', 'Number', num);
%     %     path = 'R_zero_data/Figs/';
%     %     if(~isempty(hFig))
%     %         saveas(num, [path, 'fig4', label], 'emf');
%     %         saveas(num, [path, 'fig4', label], 'eps');
%     %         saveas(num, [path, 'fig4', label], 'fig');
%     %         saveas(num, [path, 'fig4', label], 'png');
%     %     end
% end