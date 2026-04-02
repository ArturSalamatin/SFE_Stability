clc
clear all
close all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)
%% params
f = [1.87, 3.2, 5.2, 12, 40, 100];% linspace(1,505,11);
R = 0.8; %0.2;% [0.1, 0.3, 0.5,1,2,2.5];
col = {'k'};
style = {'-', '--', '-', '-.', '-'};
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
    pen = set_pen(col, style);
    out = calc_sigma(f, R, starter, pen);  
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

% figure(1)
% hold on
% y = [1.87, 3.2, 5.2, 12];
% x = ones(size(y))*2.0;
% plot(x,y, 'ok', 'MarkerFaceColor', 'black')

