clc
clear all
% close all

global DEBUG
DEBUG = true;

params.t = 0.45;
params.a = sqrt(2*params.t);
params.z2 = z2(params);


params.R = 1;
params.alpha = 5;
sigma = 100;

out = J(2.966645999, params);