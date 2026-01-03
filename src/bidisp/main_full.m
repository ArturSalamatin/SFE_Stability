clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = 0;
sigma_min_limit = -2.2;
q = 1.008;
xBarLeft = 1e-3;
xBarRight = 0*4e-2;

solver = @(problem, mesh, params) solver_KellerBox(problem, mesh, params);
label = '';
N = 200;

R_vals = linspace(0,3,31);
h_vals = [0.1, 0.5, 1, 2, 3, 5, 10];

A0_vals = 0.2;
Alpha_vals = [0.1, 0.3, 0.5, 0.7];
Tau0_vals = 0.2;

for i = 1:numel(Tau0_vals)
    a0 = A0_vals;
    tau0 = Tau0_vals(i);
    num = 5200+tau0*100;
    my_figure(num)
    
    for j = 1:numel(Alpha_vals)
        alpha = Alpha_vals(j);
        
        col = {'k'};
        style = {'-'};
        pen = set_pen(col, style);
        %% run calculations
        %     [tau0, a0, sigma] = ...
        %         calc_sigma_R_zero(solver, N, alpha, tau0, a0, pen);
        %     save(['full_data/alpha_', num2str(alpha*100), '.mat'], ...
        %         'tau0', 'a0', 'sigma', 'alpha')
    end
end


