function sol = tr_starter_R_zero_X_Psi(solver, sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
% sol = transform(solver(problem, mesh), params, problem.base_state, sigma);
sol = transform(...
    solver(problem, mesh), params, sigma, mesh);
end

function problem = set_problem(mesh, sigma, params)
problem.eqN = 2; % number of equations
problem.ids = 1:problem.eqN; % iterator for rows/cols within a block

% problem.base_state = set_base_state(mesh, params);

% c = problem.base_state.mid_c;
% xBar = problem.base_state.mid_x/params.a;
% mid_x = problem.base_state.mid_x;
% mid_xi = z_of_x(mid_x, params)/params.z2;
% figure(80)
% plot(mid_xi, (1-c)./(xBar.*xBar))
% figure(81)
% plot(mid_xi, (1-c))
% figure(82)
% plot(mid_xi, (1-c)./xBar)

problem.M = problem.eqN * mesh.N; % nmbr of discrete unknows
% problem.block_matrix = @(xL, xR, segm_i)block(...
%     xL, xR, segm_i, sigma, problem);
% problem.diag = @(i)Diag(i, problem.eqN, problem.base_state);
% problem.BC = @()BC_Psi_LR(params, problem.base_state, sigma, problem.eqN, mesh);
problem.BC_L = @()BC_L(params, [], sigma, problem.eqN, mesh);
problem.BC_R = @()BC_R(params, [], sigma, problem.eqN, mesh);
% problem.JC = @()JC(params, problem.base_state, problem.eqN);
problem.RK = @(x,y)my_ode(x,y, sigma, params);
end

function out = BC_L(params, ~, sigma, eqN, mesh)
% [Phi, Z]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0

x_left = mesh.xBarL*params.a;

[~,~, Phi_left, Z_left] = tr_calc_inlet_solution_assymptotics(...
    x_left, params, sigma);
%% BC at the left end
left = eye(eqN,eqN);
% left(1,:) = [0,1]; % X(left) = X(x_s)
% left(1,1) = 1; % Psi(left) = Psi_left, left = delta -> 0
% left(2,2) = 1; % X(left) = X_left, left = dalta -> 0
%% BC at the right end
right = zeros(eqN,eqN);
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
% X_left = Y_left/xi_left;
out.rhs = [Phi_left; Z_left];
% out.rhs = [Psi_left;X_left];
end

function out = BC_R(params, ~, sigma, eqN, ~)
% [Phi, Z]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
a = params.a;
a0 = params.a0;
g1 = params.g1;
g0 = params.g0 - g1;
dz0dt = params.dz0dt;
xi0 = params.xi0;

%% BC at the left end
left = zeros(eqN,eqN);
%% BC at the right end
right = eye(eqN,eqN);
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
Phi_right = -((g1+g0*(1-xi0))*a*dz0dt + (1+sigma)*a0/a)/(a-a0);
Z_right = 1;
out.rhs = [Phi_right, Z_right];
end

function sol = transform(sol, params, sigma, mesh)
% [Phi, Z]
global xBarRight
t = reshape(sol.t, numel(sol.t), 1);% sol.t';
a0 = params.a0;
a = params.a;
dz0dt = params.dz0dt;
g1 = params.g1;
x_right = a0*(1+xBarRight);

r = (2+sigma)/params.C2;

[x,xi, Psi_r, X_r] = calc_solution_assymptotics(...
    x_right, params, sigma);
%% transform
% get Psi from Phi
sol.y(:,1) = sol.y(:,1).*((a-mesh.x').^(r+1));
% get X from Z
sol.y(:,2) = sol.y(:,2).*((a-mesh.x').^(r));
%% normalize
ff = sol.y(end,2)/X_r;
sol.factor = ff;% sol.y(end,2)/X_r
sol.y = sol.y/sol.factor;
%% input
%[Psi, X]
Psi = sol.y(:,1);
X = sol.y(:,2);
Phi = 0*Psi;
Gamma = 0*Psi;
%% output
%[Psi, X, Phi, Gamma, Omega, Y, Psi+X]
Y = X.*t;
% Omega
g0 = params.g0;
dz2dt = params.dz2dt;
a = params.a;
Omega = Psi + Y*a*g0*dz2dt;
%% flux Q, pressure P
Q = 0*(Phi+Psi)./(1-t);
% calc P
R = params.R;
mu = exp(R*t); % viscosity
grad_p = -mu.*(Psi./t + Phi);
grad_p(1)=0;
I = 1:numel(mu)-1;
P = (grad_p(I)+grad_p(I+1))/2.*(t(I+1)-t(I));
P = [0; cumsum(P)];
dP = grad_p;
%% form return variable
sol.y = [Psi, X, Phi, Gamma, Omega, Y, 0*(Psi+X), 0*Q, 0*P, 0*dP];
%% Keller-box scheme behind the jump point
C2 = params.C2;
C1 = params.C1;
% r = 1+(2+sigma)/C2;
a = params.a;
g1 = params.g1;
g0 = params.g0 - g1;
dz0dt = params.dz0dt;
xi0 = params.xi0;
% sol.condition = Psi(1); % Psi(xi=0) == 0
sol.condition = ...
    -(g1+g0*(1-xi0))*a*dz0dt - (1+sigma)*a0/a; % Psi(xi=xi0) == -(g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0)))

% x_right = a0*(1+a0/5);
% [x,xi, Psi, X] = calc_solution_assymptotics(...
%     x_right, params, sigma);
sol.condition = Psi_r;
%% RK integration to the left from the jump point
% C2 = base_state.C2;
% C1 = base_state.C1;
% a = params.a;
% dz0dt = base_state.dz0dt;
% g0 = params.g0 - params.g1;
% g1 = params.g1;
% xi0 = base_state.z0/base_state.z2;
% sol.BC = Psi(end) - X(end)*(g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0)));
%% RK integration FROM the jump point to the left
% a = params.a; % == sqrt(2tau)
% tau = (a*a)/2.0;
% alpha = params.r;
% C2 = base_state.C2;
% C1 = base_state.C1;
% zeta2 = base_state.z2;
% xi_left = mesh.left.L;
% [Psi_left, X_left] = ...
%     calculate_X_Psi(xi_left, tau, alpha, C1, C2, zeta2, sigma);
% sol.BC = Psi_left/X_left - Psi(1)/X(1);
% sol.base_state = base_state;
sol.sigma = sigma;
end

function dy = my_ode(x,y, sigma, params)
% [Phi, Z]
alpha = params.r;
a = params.a; % == sqrt(2tau)
g1 = 1-alpha;
Da =  alpha+(1-alpha)*a;
Dx = alpha+(1-alpha)*x;
z = z_of_x(x, params);
C2 = params.C2;

r = (2+sigma)/C2;

f = zeros(2,2);
f(1,1) = g1/Dx +(r+1)/(a-x);
f(1,2) = g1*a/(x*Da*(a-x));

f(2,1) = -a*(a-x)/(Dx*z*C2);
f(2,2) = -(a*a/(x*Da) + x*(1+sigma)/Dx - C2*r*z/(a-x))/(z*C2);

dy = f*y;
end


