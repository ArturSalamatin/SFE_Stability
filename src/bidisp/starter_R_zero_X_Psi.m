function sol = starter_R_zero_X_Psi(solver, sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
% sol = transform(solver(problem, mesh), params, problem.base_state, sigma);
sol = transform(...
    solver(problem, mesh), params, sigma, [] ...problem.base_state
    );
end

function problem = set_problem(mesh, sigma, params)
problem.eqN = 2; % number of equations
problem.ids = 1:problem.eqN; % iterator for rows/cols within a block

% problem.base_state = set_base_state(mesh, params);

problem.M = problem.eqN * mesh.N; % nmbr of discrete unknows
% problem.block_matrix = @(xL, xR, segm_i)block(...
%     xL, xR, segm_i, sigma, problem);
% problem.diag = @(i)Diag(i, problem.eqN, problem.base_state);
% problem.BC = @()BC_Psi_LR(params, problem.base_state, sigma, problem.eqN, mesh);
problem.BC_L = @()BC_L(params, [] ...problem.base_state
    , sigma, problem.eqN, mesh);
problem.BC_R = @()BC_R(params, [] ...problem.base_state
    , sigma, problem.eqN, mesh);
problem.JC = @()JC(params, [] ...problem.base_state
    , problem.eqN);
problem.RK = @(x,y)my_ode(x,y, sigma, params);

problem.BVP_f = @(u,y) sing_bvp_ode(u,y,sigma,params);
problem.BVP_S = @() S(sigma,params.C2);
problem.BVP_bc = @(ya,yb) bvp_bc_fcn(ya,yb, mesh, sigma, params);

problem.sigma = sigma;
end

function base_state = set_base_state(mesh, params)
a = params.a;

I = 1:(mesh.N-1);
nodes_xBar = mesh.xBar; % mesh nodes
mid_nodes_xBar = (nodes_xBar(I)+nodes_xBar(I+1))/2; % centers of mesh segments

base_state.params = params;

base_state.mid_x = mid_nodes_xBar*a;
base_state.mid_c = c_of_x(base_state.mid_x, params);
base_state.mid_g = g(base_state.mid_x, params);
base_state.mid_dcdz = dcdz(base_state.mid_x, params);
base_state.mid_dxBardt = dxBardt(base_state.mid_x, params);
end

function out = block(xiL, xiR, segm_i, sigma, problem)
% [Psi, X]
eqN = problem.eqN;
a = problem.base_state.params.a;
%% set params
% def: C1(t) == sqrt(2t)/zeta2
C1 = problem.base_state.params.C1;
% def: C2(t) == 2t*dzeta2dt/zeta2
C2 = problem.base_state.params.C2;

c = problem.base_state.mid_c(segm_i);
g = problem.base_state.mid_g(segm_i);
xBar = problem.base_state.mid_x(segm_i)/a;

% segm_nmbr = numel(segm_i);

xiMid  =(xiL+xiR)/2;
%% set out
out = zeros(eqN,eqN);
out(1,1) = -g/C1./xBar;
out(1,2) = -g/C1*(1-c)./(xBar.^2);

out(2,1) = 1./(C2.*xBar.*xiMid);
out(2,2) = ((1-c)./(xBar.^2)+1+sigma)./(xiMid*C2);
end

function out = Diag(~, eqN, ~)
out = sparse(1:eqN, 1:eqN, ones(1,eqN), eqN, eqN, eqN);
end

function out = BC_L(params, base_state, sigma, eqN, mesh)
% [Psi, X]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
C2 = params.C2;
C1 = params.C1;
% r = 1+(2+sigma)/C2;
a = params.a;
g1 = params.g1;
g0 = params.g0 - g1;
dz0dt = params.dz0dt;
xi0 = params.xi0;

x_left = mesh.xBarL*params.a;

% [~, X_left] = ...
%     calculate_X_Psi(xi_left, tau, alpha, C1, C2, zeta2, sigma);
[~,~, Psi_left, X_left] = inlet_solution_asymptotics(...
    x_left, params, sigma);
%% BC at the left end
left = eye(eqN,eqN);
% left(1,:) = [0,1]; % X(left) = X(x_s)
% left(1,1) = 1; % Psi(left) = Psi_left, left = delta -> 0
% left(2,2) = 1; % X(left) = X_left, left = dalta -> 0
%% BC at the right end
right = zeros(eqN,eqN);
% Psi(end) - X(end)*(g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0))) = 0
% right(2,:) = [1, -(g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0)))];
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
% X_left = Y_left/xi_left;
out.rhs = [Psi_left; X_left];
% out.rhs = [Psi_left;X_left];
end

function out = BC_Psi_LR(params, base_state, sigma, eqN, mesh)
% [Psi, X]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
a = params.a;
a0 = params.a0;
g1 = params.g1;
dz0dt = params.dz0dt;

x_left = mesh.xBarL*params.a;
[~,~, Psi_left, ~] = inlet_solution_asymptotics(...
    x_left, params, sigma);


left = zeros(eqN,eqN);
right = zeros(eqN,eqN);
rhs = zeros(eqN,1);
%% BC at the left end
% 1*Psi(left) +0*Psi(right) = 0
left(1,:) = [1,0];
right(1,:) = [0,0];
rhs(1) = Psi_left;
%% BC at the right end
% 0*Psi(left) + 1*Psi(right) = Psi_r == = -g1*a*dz0dt - (1+sigma)*a0/a
left(2,:) = [0,0];
right(2,:) = [1,g1*a*dz0dt + (1+sigma)*a0/a];
rhs(2) = 0;

out.left = left;
out.right = right;
out.rhs = rhs;
end

function out = BC_R(params, base_state, sigma, eqN, mesh)
% [Psi, X]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
C2 = params.C2;
C1 = params.C1;
% r = 1+(2+sigma)/C2;
a = params.a;
a0 = params.a0;
g1 = params.g1;
dz0dt = params.dz0dt;

%% BC at the left end
left = zeros(eqN,eqN);
% left(1,:) = [0,1]; % X(left) = X(x_s)
% left(1,1) = 1; % Psi(left) = Psi_left, left = delta -> 0
% left(2,2) = 1; % X(left) = X_left, left = dalta -> 0
%% BC at the right end
right = eye(eqN,eqN);
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
out.rhs = [-g1*a*dz0dt - (1+sigma)*a0/a, 1];
end

function out = JC(params, base_state, eqN)
%[Psi, X]

% left*y(left) + right*y(right) = rhs
% [Psi] + a*dz0dt*g0*X = 0, Y = xi0*X
% - 1*Psi(left) + 1*Psi(right) + a*dz0dt*g0*X(right) = 0
% - 1*X(left)   + 1*X(right) = 0
a = params.a; % a == sqrt(2*t)
g0 = params.g0-params.g1;
dz2dt = params.dz2dt;
out.left = -eye(eqN,eqN);
out.right = eye(eqN,eqN);
out.right(1,2) = a*dz2dt*g0;
out.rhs = [0;0];
end

function sol = transform(sol, params, sigma, base_state)
% global xBarRight
t = reshape(sol.t, numel(sol.t), 1);% sol.t';
a0 = params.a0;
a = params.a;
dz0dt = params.dz0dt;
g1 = params.g1;
x_right = a0;

[x,xi, Psi_r, X_r] = calc_solution_assymptotics(...
    x_right, params, sigma);

%% normalize
ff = sol.y(end,2)/X_r;
sol.factor = 1;% ff;
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
    -(g1*a*dz0dt + (1+sigma)*a0/a)*X(sol.id); % Psi(xi=xi0) == -(g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0)))

% x_right = a0*(1+a0/5);
% [x,xi, Psi, X] = calc_solution_assymptotics(...
%     x_right, params, sigma);
% sol.condition = Psi_r;
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
% [Psi, X]
a = params.a; % == sqrt(2tau)
z = z_of_x(x, params);
C2 = params.C2;

gx = g(x, params);
Gx = G_of_x(x, params);
Ga = G_of_x(a, params);

f = zeros(2,2);
f(1,1) = gx/Gx;
f(1,2) = gx*a/(x*Ga);

f(2,1) = -a/(Gx*z*C2);
f(2,2) = -(a*a/(x*Ga) + x*(1+sigma)/Gx)/(z*C2);

dy = f*y;
end

function out = S(sigma, C2)

out = [0, 0; -[1, 2+sigma]/C2];

end

function dy = sing_bvp_ode(u,y, sigma, params)
% [Psi, X]
dy = -my_ode(params.a - u,y, sigma, params);
end

function out = bvp_bc_fcn(ya,yb, mesh, sigma, params)
global xBarLeft

x_left = mesh.x(1);

[~,~, Psi, X] = inlet_solution_asymptotics(...
    x_left, params, sigma);

out = [ya(1)*X-ya(2)*Psi, yb(2)-1];
% out = [ya(1) - Psi, ya(2) - X];
end


