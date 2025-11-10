clc
clear all
% close all

global sigma_fig
sigma_fig = 6;

%% packed bed params
a0 = 0.2;
alpha = [0, 0.1, 0.3, 0.5, 0.7];
marker = ['o','s','d','*','x'];

params = poly_case(a0, alpha(4), 0.1);
params.marker = marker(2);

N = 900;
q = 1.06;
mesh = log_mesh(q, N, params);
starter = @(sigma, params) starter_Y(...
                sigma, params, mesh);

col = {'k'};
style = {'-'};
pen = set_pen(col, style);

h = 10;
R = 0.1;
[out, sol] = calc_sigma(h, R, starter, pen, params);
