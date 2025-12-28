clc
% close all
clear all

[params, sigma] = case_2();
N = 7001;
mesh = quasiuniform_mesh(N, params);
% mesh = log_mesh(1.0005, N, params);

sol = starter_Y(sigma, params, mesh);

pen.lc = 'b';
pen.style = '-';
plot_solution(sol.t,sol.y,pen);



function [params, sigma] = case_1
params.a0 = 0.3;
params.a1 = 1.0;
params.r = 0.5; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.0001;
params.h = 5;
params.t = 0.4;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);

sigma = -1.62442043;
end

function [params, sigma] = case_2
params.a0 = 0.01;
params.a1 = 1.0;
params.r = 0.2; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.0001;
params.h = 5;
params.t = 0.47;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);

sigma = -2.061672191741087;
end
