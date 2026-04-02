clc
clear all
close all

%% params
h = [15, 0.5, 1, 2, 3, 5, 10];% linspace(1,505,11);
R = 0.1:0.1:1.5;
R  = R(1);
h = h(1);
col = {'m'};
style = {'-', '--', '-', '-.', '-'};
[params, ~] = case_monodisp;
%% set mesh
N = [2001];
for i = 1:numel(N)
    %% choose starter
%     starter = @(sigma, params) starter_omega_X(sigma, params, mesh);
%     pen = set_pen(col, style{2});
%     out = calc_sigma(f, R, starter, pen);
    %% choose starter
    starter = @(sigma, params) starter_Y(sigma, params, N(i));
    pen = set_pen(col, style);
    out = calc_sigma(h, R, starter, pen, params);  
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

function [params, sigma] = case_1
params.a0 = 0.1; %0.5;
params.a1 = 1.0;
params.r = 0.05; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.5;
params.h = 5;
params.t = 0.2;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);
params.dz2dt = dz2dt(params);
params.C1 = params.a./params.z2;
params.C2 = params.C1*params.a.*params.dz2dt;

sigma = -1.9072018;
end

function [params, sigma] = case_monodisp
params.a0 = 0.5;
params.a1 = 1.0;
params.r = 0.0; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.5;
params.h = 5;
params.t = 0.2;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);
params.dz2dt = dz2dt(params);
params.C1 = params.a./z2(params);
params.C2 = params.C1*params.a.*params.dz2dt;

sigma = -1.9072018;
end
