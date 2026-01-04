clc
clear all
close all

global sigma_fig sigma_max_limit sigma_min_limit q xBarLeft xBarRight
sigma_fig = 9;
sigma_max_limit = 5;
sigma_min_limit = -2.2;
q = 1.008;
xBarLeft = 1e-3;
xBarRight = 0*4e-2;

solver = @(problem, mesh, params) solver_KellerBox(problem, mesh, params);
label = '';
N = 200;

R_vals = [0.1, 0.3, 1, 2];
h_vals = linspace(0.01, 40, 101);

A0_vals = 0.2;
Alpha_vals = [0.1, 0.3, 0.5, 0.7, 0.9];
Tau0_vals = [0.03, 0.2];

% cols = {'m', 'k', 'b', 'r', [0.1 0.1 0.1]};

cols = {[0, 0.4470, 0.7410]
    'red'
    [0.9290, 0.6940, 0.1250]
    'blue'
    [0.4660, 0.6740, 0.1880]
    [0.3010, 0.7450, 0.9330]
    [0.6350, 0.0780, 0.1840]};

for i = 1:numel(R_vals)
    R = R_vals(i);
    a0 = A0_vals;
    tau0 = Tau0_vals;
    num = 5000+R*100;
    my_figure(num)
    hold on
    axis([0 max(h_vals) -2 0.5])
    set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
    set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])
    ylabel(['{\it\sigma}, ' char(8211)])
    xlabel(['{\ith}, ' char(8211)])
    
    Sigma = zeros(numel(Alpha_vals), numel(Tau0_vals), numel(h_vals));
        load(['full_data/R_', num2str(R*100), '.mat'], ...
            'R_vals', 'Sigma');
    tic
    for j = 1:numel(Alpha_vals)
        alpha = Alpha_vals(j);
        disp(['alpha = ', num2str(alpha)]);
        col = {'k'};
        style = {'-'};
        pens = set_pen(col, style);
        %% run calculations
%         Sigma(j,:,:) = ...
%             calc_sigma_full_2(solver, N, alpha, tau0, a0, R, h_vals, pens);
        sigma = reshape(Sigma(j,:,:), size(Sigma, [2,3]));
        
        plot(h_vals, sigma(1,:), '-',  'Color', cols{j}, 'LineWidth', 1)
        plot(h_vals, sigma(2,:), '--', 'Color', cols{j}, 'LineWidth', 1)
    end
    toc
%     save(['full_data/R_', num2str(R*100), '.mat'], ...
%         'R_vals', 'Sigma')
    
    label = '';
    switch R
        case 0.1
            label = 'A';
            annotation('textbox',...
                [0.178795331790121 0.782747603833865 0.0868296682098794 0.115015974440895],...
                'String',{label},...
                'LineStyle','none',...
                'FontSize',16,...
                'FontName','Times New Roman',...
                'FitBoxToText','off',...
                'BackgroundColor',[1 1 1]);
        case 0.3
            label = 'B';
            annotation('textbox',...
                [0.178795331790121 0.782747603833865 0.0868296682098794 0.115015974440895],...
                'String',{label},...
                'LineStyle','none',...
                'FontSize',16,...
                'FontName','Times New Roman',...
                'FitBoxToText','off',...
                'BackgroundColor',[1 1 1]);
        case 1
            label = 'C';
            annotation('textbox',...
                [0.81942033179012 0.233226837060702 0.0868296682098796 0.115015974440895],...
                'String',{label},...
                'LineStyle','none',...
                'FontSize',16,...
                'FontName','Times New Roman',...
                'FitBoxToText','off',...
                'BackgroundColor',[1 1 1]);
        case 2
            label = 'D';
            annotation('textbox',...
                [0.81942033179012 0.233226837060702 0.0868296682098796 0.115015974440895],...
                'String',{label},...
                'LineStyle','none',...
                'FontSize',16,...
                'FontName','Times New Roman',...
                'FitBoxToText','off',...
                'BackgroundColor',[1 1 1]);
    end
    
%     hFig = findobj('Type', 'figure', 'Number', num);
%     path = 'full_data/Figs/';
%     figid = 'fig8';
%     if(~isempty(hFig))
%         saveas(num, [path, figid, label], 'emf');
%         saveas(num, [path, figid, label], 'eps');
%         saveas(num, [path, figid, label], 'fig');
%         saveas(num, [path, figid, label], 'png');
%     end
end


