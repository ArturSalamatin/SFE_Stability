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
N = 1000;

R_vals = linspace(0,3,101);
h_vals = [0.1, 0.5, 1, 2, 3, 5, 10];

A0_vals = 0.2;
Alpha_vals = 0.1; % [0.1, 0.3, 0.5, 0.7];
Tau0_vals = 0.2;

for i = 1:numel(Tau0_vals)
    a0 = A0_vals;
    tau0 = Tau0_vals(i);
    num = 5200+tau0*100;
    my_figure(num)
    hold on
    set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
    set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])
    ylabel(['{\it\sigma}, ' char(8211)])
    xlabel(['{\itR}, ' char(8211)])
    
    Sigma = zeros(numel(Alpha_vals), numel(R_vals), numel(h_vals));
%     load(['full_data/tau_', num2str(tau0*100), '.mat'], ...
%         'tau0', 'A0_vals', 'Sigma', 'Alpha_vals', 'R_vals', 'h_vals')
    for j = 1:numel(Alpha_vals)
        alpha = Alpha_vals(j);
        
        col = {'k'};
        style = {'-'};
        pen = set_pen(col, style);
        %% run calculations
        Sigma(j,:,:) = ...
            calc_sigma_full(solver, N, alpha, tau0, a0, R_vals, h_vals, pen);
        
        plot(R_vals, sigma, 'k--')
        
    end
    save(['full_data/tau_', num2str(tau0*100), '.mat'], ...
        'tau0', 'A0_vals', 'Sigma', 'Alpha_vals', 'R_vals', 'h_vals')
end


