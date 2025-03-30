clc
close all
clear all



params.a0 = 0.2;
params.a1 = 1.0;
params.r = 0.5; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.1;
params.h = 0.5;
params.t = 0.2;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);

N = 301;
mesh = quasiuniform_mesh(1, params.a0, 0, N, params);

sigma = -1.9;
sol = starter_Y(sigma, params, mesh);

pen.lc = 'black';
pen.style = '-';
plot_solution(sol.t,sol.y,pen);