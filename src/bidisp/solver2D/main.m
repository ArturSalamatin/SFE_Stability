clc
close all
clear all

% scripts and functions from other folder can be run
% addpath('../')

global x_tol
x_tol = 1e-10;
%% set parameters
H = 1;
B = 0;
R = 0.0;
a0 = 0.1;
alpha = 0.5;
Nr = 2;
Nz = 365;
dt = 0.01;
T = 0.5;
params = set_params(a0, alpha, R, B, H);
%% set mesh
mesh = set_mesh(Nr, Nz, dt, T, params);

states = solver2D(mesh, params);

L = 1:(mesh.size)/2;
p = states(0*mesh.size+L,:);
c = states(1*mesh.size+L,:);
y = states(2*mesh.size+L,:);
G = states(3*mesh.size+L,:);
x = states(4*mesh.size+L,:);

if(max(max(abs(y - x.*x/2))) > 1e-10)
    error("y is incorrect!");
end
if(max(max(abs(G_of_x(x, params) - G))) > 1e-10)
    error("G is incorrect!");
end


figure(1001)
hold on
plot(mesh.z, x(:, end), 'k')
axis([0 H 0 1])

z = z_of_x(x(:,end), params, T);
plot(z, x(:, end), 'r')
axis([0 H 0 1])

figure(1002)
hold on
plot(mesh.z, c(:, end), 'k')
axis([0 H 0 1])

z = z_of_x(x(:,end), params, T);
C = c_of_x(x(:,end), params, T);
plot(z, C, 'r')
axis([0 H 0 1])


