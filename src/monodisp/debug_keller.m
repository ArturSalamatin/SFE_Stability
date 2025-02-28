clc
clear all
% close all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)
%% params
f = 5; %[0.1,0.5,1,5,10];% linspace(1,505,11);
R = 0.2;% [0.1, 0.3, 0.5,1,2,2.5];
col = {'b'};
style = {'-', '--'};
%% set mesh
Left = 0;
Right = 1;
N = [5001];
for i = 1:numel(N)
    mesh = uniform_mesh(Left, Right, N(i));
    %% choose starter
%     starter = @(sigma, params) starter_omega_X(sigma, params, mesh);
%     pen = set_pen(col, style{2});
%     out = calc_sigma(f, R, starter, pen);
    %% choose starter
    starter = @(sigma, params) starter_omega_Y(sigma, params, mesh);
    pen = set_pen(col, style{1});
    out = calc_sigma(f, R, starter, pen);  
    %% choose starter
%     A = 35;
%     mesh = uniform_mesh(-A, A, N(i));
%     starter = @(sigma, params) starter_inf_omega_Y(sigma, params, mesh);
%     pen = set_pen(col, style{2});
%     out = calc_sigma(f, R, starter, pen);  
end


figure(8)
hold on
box on
xlabel('{\itR}')
ylabel('{\it\sigma}')
plot(R, out, 'k-', 'LineWidth', 1)



