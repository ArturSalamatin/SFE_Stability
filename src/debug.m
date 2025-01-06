clc
clear all
close all

params.a0 = 0.05;
params.r = 0.5;
params.t = (params.a0^2)/2;
params.a = sqrt(2*params.t);
params.g0 = 1-params.r + params.r/params.a0;
params.g1 = 1-params.r;
params.z0 = z0(params);
params.z2 = z2(params);

x = x_grid(params);
z = z_of_x(x, params);

plot(z/params.z0,x/params.a)



