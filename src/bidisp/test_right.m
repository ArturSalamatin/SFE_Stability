clc
close all

global q
q = 1.00001;

alpha = 0.2;
a0 = 0.1;
tau0 = 0.2864;
sigma = -1.63;

params = poly_case(a0, alpha, tau0, R);
g1 = params.g1;
g0 = params.g0 - g1;
xi0 = params.z0/params.z2;
C1 = params.C1;
C2 = params.C2;
psi0 = -(g0+g1)/C1*(C2+(1+sigma)*(1-xi0));
y0 = [psi0;1];


N = 250;
delta = a0/10;
mesh = set_right_mesh(N, params, delta);


solver = @(problem, mesh) solver_RK(problem, mesh, params);
% starter = @(sigma, params) starter_R_zero_X_Psi(...
%     solver, sigma, params, mesh);
starter = @(sigma, params) starter_R_zero_X_Psi(...
    solver, sigma, params, mesh);

sol = starter(sigma, params);
pen = set_pen('k', '--');
plot_solution(pen, params, sol)
