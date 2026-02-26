clc
close all
clear all


% scripts and functions from other folder can be run
% addpath('../')

%% set parameters
H = 10;
B = 50;
a0 = 0.1;
alpha = 0.5;
Nr = 20;
Nz = 50;
params = set_params(a0, alpha, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, params);
%% set initial conditions
IC = initial_conditions(mesh, params);