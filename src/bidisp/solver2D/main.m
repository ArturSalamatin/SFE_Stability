clc
close all
clear all

% scripts and functions from other folder can be run
% addpath('../')

global x_tol
x_tol = 1e-10;
%% set parameters
B = 3;
R = 1.0;
a0 = 0.2;
alpha = 0.7;
Nr = 90;
Nz = 150;
H = 1.5;
dt = 0.01;
T = 1.5;
params = set_params(a0, alpha, R, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, dt, T, params);

states = solver2D(mesh, params);
plotter_2D(states, mesh);
plotter_1D(states, mesh, params);
