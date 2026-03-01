clc
close all
clear all

% scripts and functions from other folder can be run
% addpath('../')

%% set parameters
H = 10;
B = 0;
R = 0.0;
a0 = 0.1;
alpha = 0.5;
Nr = 3;
Nz = 100;
dt = 0.01;
T = 0.5;
params = set_params(a0, alpha, R, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, dt, T, params);

states = solver2D(mesh, params);

L = 1:mesh.size;
c = states(mesh.size+L,:);
x = states(4*mesh.size+L,:);