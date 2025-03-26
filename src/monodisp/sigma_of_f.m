clc
clear all
% close all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)
%% params
f = [linspace(0.01, 1, 30), linspace(1, 40, 155)];% linspace(1,505,11);
R = 0.2:0.2:2.0;
col = {'b'};
style = {'-', '--'};
%% set mesh
Left = 0;
Right = 1;
N = [401];
for i = 1:numel(N)
    mesh = uniform_mesh(Left, Right, N(i));
    %% choose starter
    %     starter = @(sigma, params) starter_omega_X(sigma, params, mesh);
    %     pen = set_pen(col, style{2});
    %     out = calc_sigma(f, R, starter, pen);
    %% choose starter
    starter = @(sigma, params) starter_omega_Y(sigma, params, mesh);
    pen = set_pen(col, style{1});
    out = calc_sigma2(f, R, starter, pen);
    
    %     figure(10)
    %     [f1, R1] = meshgrid([0,R],[0,f]);
    %     out1 = zeros(size(f1))-2;
    %     out1(2:end, 2:end) = out;
    %
    %     contour(f1,R1,out1,-2:0.2:0,'ShowText','off')
    %     ylabel('{\itf}')
    %     xlabel('{\itR}')
    % annotation('arrow',[0.898571428571429 0.208571428571429],...
    %     [0.24952380952381 0.220952380952381]);
    
    
    %% choose starter
    %     A = 35;
    %     mesh = uniform_mesh(-A, A, N(i));
    %     starter = @(sigma, params) starter_inf_omega_Y(sigma, params, mesh);
    %     pen = set_pen(col, style{2});
    %     out = calc_sigma(f, R, starter, pen);
end


% figure(8)
% hold on
% box on
% xlabel('{\itR}')
% ylabel('{\it\sigma}')
% plot(R, out, 'k-', 'LineWidth', 1)



function out = calc_sigma2(f, R, starter, pens)
out = zeros(numel(f), numel(R)) - 1.999999;

for i = 1:numel(R)
    params.R = R(i);
    for j = 1:numel(f)
        params.f = f(j);
        %         I = max(1,i-1);
        J = max(1,j-1);
        
        
        pen = set_pen(pens.lc{min(i, numel(pens.lc))}, pens.style);
        params.pen = pen;
        
        [sigma] = fit_sigma(starter, params, out(J,i));
        out(j,i) = sigma;
    end
    
    %     figure(7)
    %     hold on
    %     box on
    %     xlabel('{\itf}')
    %     ylabel('{\it\sigma}')
    %     plot(f, out(:,i), 'k-', 'LineWidth', 1)
    %
    %
    figure(9)
    hold on
    box on
    xlabel('{\itf}')
    ylabel('({\it\sigma}+2)/({\it\sigma}_{{\itf}\rightarrow\infty}+2)')
    plot(f, (out(:,i)+2)/(out(end,i)+2), 'k-', 'LineWidth', 1)
    %         save_figs(sigma, params);
    
end

figure(7)
annotation('arrow',[0.8375 0.840357142857143],...
    [0.13947619047619 0.650952380952381]);

figure(9)
axis([0 20 0 1.5])
annotation('arrow',[0.3375 0.1875],[0.570952380952381 0.87]);
end

