clc
clear all
% close all

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
for alpha = [0.5, 0.2:0.2:0.8]
    num = 2200+alpha*10;
    clc
    col = {'k'};
    style = {'-'};
    pen = set_pen(col, style);
    %% make calculation grid
    eps = 1e-2;
    tau0 = linspace(eps, 0.5 - eps, 300);
    a0 = linspace(eps,1-eps, 300);
    [tau0, a0] = meshgrid(tau0, a0);
    %% run calculations
    [tau0, a0, sigma] = ...
        calc_sigma_R_zero(solver, N, alpha, tau0, a0, pen);
    save(['R_zero_data/alpha_', num2str(alpha*100), '.mat'], ...
        'tau0', 'a0', 'sigma', 'alpha')
    
        figure(num)
        hold on
        axis([0 0.5 0 1])
        contour(tau0, a0, sigma, ...
            ...linspace(-2,-1,11),...
            'ShowText','on',...
            'linecolor', 'black')
    
        x = linspace(0, 0.5, 700);
        y = sqrt(2*x);
        plot(x,y, 'k-', 'linewidth', 1)
    
        x = -1:0.1:0;
        y = -(1+x)*alpha./(x*(1-alpha));%1/(1-alpha)*(alpha./(x+2)-alpha);
        z = y.*y/2;
        plot(z,y,'ks', 'MarkerFaceColor', 'black')
    
        switch alpha
            case 0.2
                label = 'A';
            case 0.4
                label = 'B';
            case 0.6
                label = 'C';
            case 0.8
                label = 'D';
        end
    
        box on
        set(gca, 'xTick', [0:0.1:0.5])
        set(gca, 'yTick', [0:0.2:1.0])
        set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
        set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])
        ylabel(['{\ita}_0, ' char(8211)])
        xlabel(['{\it\tau}_0, ' char(8211)])
        annotation('textbox',...
            [0.207441165123454 0.793200250350447 0.115740738211223 0.141129029253798],...
            'String',{label},...
            'LineStyle','none',...
            'FontSize',16,...
            'FontName','Times New Roman');
    
        hFig = findobj('Type', 'figure', 'Number', num);
        path = 'R_zero_data/Figs/';
        if(~isempty(hFig))
            saveas(num, [path, 'fig4', label], 'emf');
            saveas(num, [path, 'fig4', label], 'eps');
            saveas(num, [path, 'fig4', label], 'fig');
            saveas(num, [path, 'fig4', label], 'png');
        end
end
            
            
