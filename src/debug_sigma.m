clc
clear all
% close all

global DEBUG
DEBUG = true;

params.a0 = 0.4;
params.r = 0.5;
params.g1 = 1-params.r;
params.g0 = params.g1 + params.r/params.a0;

params.t = 0.45;
params.a = sqrt(2*params.t);

params.x = x_grid(params, 1e-3);
params.z = z_of_x(params.x, params);
params.c = c_of_x(params.x, params);

params.z0 = z0(params);
params.z2 = z2(params);


params.R = 1;
params.alpha = 5;
sigma = 100;

out = J(2.966645999, params);