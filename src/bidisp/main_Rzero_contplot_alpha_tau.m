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

solver = @(problem, mesh, params) solver_RK(problem, mesh, params);        
label = '';
N = 200;
for A0 = [3]
    a0 = A0/10;
    num = 2200+a0*10;
    clc
    col = {'k'};
    style = {'-'};
    pen = set_pen(col, style);
    %% make calculation grid
    eps = 1e-2;
    tau0 = linspace(a0*a0/2 + eps, 0.5 - eps, 10);
    alpha = linspace(eps,1-eps, 10);
    [tau0, alpha] = meshgrid(tau0, alpha);
    sigma = 0*tau0 -2;
    %% run calculations
%     [tau0, a0, sigma] = ...
%         calc_sigma_R_zero(solver, N, alpha, tau0, a0, pen);
%     save(['R_zero_data/a0_', num2str(a0*100), '.mat'], ...
%         'tau0', 'a0', 'sigma', 'alpha')
% load(['R_zero_data/alpha_', num2str(alpha*100), '.mat'], ...
%         'tau0', 'a0', 'sigma', 'alpha')

%     levels = linspace(-2,-1,21);
%     switch Alpha
%             case 2
%             levels = linspace(-2,-1,11);
%     end
        figure(num)
        hold on
        axis([floor((a0*a0/2)*10)/10 0.5 0 1])
        contour(tau0, alpha, sigma, ...
            ...levels,...
            'ShowText','on',...
            'linecolor', 'black')
    
        x = [1,1]*a0*a0/2;
        y = [0,1];
        plot(x,y, 'k-', 'linewidth', 1)
    
%         x = -1:0.1:0;
%         y = -(1+x)*alpha./(x*(1-alpha));%1/(1-alpha)*(alpha./(x+2)-alpha);
%         z = y.*y/2;
%         plot(z,y,'ks', 'MarkerFaceColor', 'black')
    
        label = '';
        switch A0
            case 1
                label = 'A';
            case 2
                label = 'B';
            case 3
                label = 'C';
            case 4
                label = 'D';
        end
    
        box on
%         set(gca, 'xTick', [0:0.1:0.5])
        set(gca, 'yTick', [0:0.2:1.0])
        set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
        set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])
        ylabel(['{\it\alpha}, ' char(8211)])
        xlabel(['{\it\tau}_0, ' char(8211)])
        annotation('textbox',...
            [0.845461998456785 0.809174691245016 0.115740738211223 0.141129029253798],...
            'String',{label},...
            'LineStyle','none',...
            'FontSize',16,...
            'FontName','Times New Roman');
    
        hFig = findobj('Type', 'figure', 'Number', num);
        path = 'R_zero_data/Figs/';
        if(~isempty(hFig))
            saveas(num, [path, 'fig5', label], 'emf');
            saveas(num, [path, 'fig5', label], 'eps');
            saveas(num, [path, 'fig5', label], 'fig');
            saveas(num, [path, 'fig5', label], 'png');
        end
end
            
            
