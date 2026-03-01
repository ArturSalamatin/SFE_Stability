clc
close all
clear all

% scripts and functions from other folder can be run
% addpath('../')

%% set parameters
H = 10;
B = 50;
R = 0.1;
a0 = 0.1;
alpha = 0.5;
Nr = 20;
Nz = 50;
tau = 0.01;
params = set_params(a0, alpha, R, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, params);
L = linear_index(mesh.J, mesh.I, mesh);
%% set initial conditions
IC = initial_conditions(mesh, params);
state = make_state(IC);

state_s = single_iteration(state, state, mesh, params, tau);


