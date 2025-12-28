clc
close all
clear all

[params, sigma] = case_1();
N = 3001;
mesh = quasiuniform_mesh(N, params);

sol = starter_Y(sigma, params, mesh);

pen.lc = 'black';
pen.style = '-';
plot_solution(sol.t,sol.y,pen);


function [params, sigma] = case_1
params.a0 = 0.5;
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

sigma = -1.9072018;
end