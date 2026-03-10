clc
close all
clear all

% scripts and functions from other folder can be run
% addpath('../')

global x_tol
x_tol = 1e-10;
%% set parameters
H = 1.5;
B = 200;
R = 1.8;
a0 = 0.1;
alpha = 0.3;
Nr = 10;
Nz = 50;
dt = 0.002;
T = 0.5;
params = set_params(a0, alpha, R, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, dt, T, params);

states = solver2D(mesh, params);
plotter_2D(states, mesh);
plotter_1D(states, mesh, params);
