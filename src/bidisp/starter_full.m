function sol = starter_full(solver, sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
% sol = transform(solver(problem, mesh), params, problem.base_state, sigma);
sol = transform(...
    solver(problem, mesh), params, sigma, [] ...problem.base_state
    );
end

function problem = set_problem(mesh, sigma, params)
problem.eqN = 4; % number of equations
problem.ids = 1:problem.eqN; % iterator for rows/cols within a block

problem.base_state = set_base_state(mesh, params);

problem.M = problem.eqN * mesh.N; % nmbr of discrete unknows
problem.block_matrix = @(segm_i)block_matrix(...
    segm_i, sigma, problem.base_state, problem.eqN);
problem.diag = @(i)Diag(i, problem.eqN, problem.base_state);
problem.BC_Keller = @()BC_Keller(params, problem.eqN);
problem.BC_L = @()BC_L(params, [] ...problem.base_state
    , sigma, problem.eqN, mesh);
% problem.BC_R = @()BC_R(params, [] ...problem.base_state
%     , sigma, problem.eqN, mesh);
problem.JC = @()JC(params, [] ...problem.base_state
    , problem.eqN);
problem.RK = @(x,y)my_ode(x,y, sigma, params);

% problem.BVP_f = @(u,y) sing_bvp_ode(u,y,sigma,params);
% problem.BVP_S = @() S(sigma,params.C2);
% problem.BVP_bc = @(ya,yb) bvp_bc_fcn(ya,yb, mesh, sigma, params);

problem.sigma = sigma;
end

function base_state = set_base_state(mesh, params)
a = params.a;

I = 1:(mesh.N-1);
nodes_xBar = mesh.xBar; % mesh nodes
nodes_xi = mesh.xi; % mesh nodes
mid_nodes_xBar = (nodes_xBar(I)+nodes_xBar(I+1))/2; % centers of mesh segments
mid_nodes_xi = (nodes_xi(I)+nodes_xi(I+1))/2; % centers of mesh segments

base_state.params = params;

base_state.mid_xBar = mid_nodes_xBar;
base_state.mid_x = mid_nodes_xBar*a;
base_state.mid_xi = mid_nodes_xi;
base_state.mid_c = c_of_x(base_state.mid_x, params);
base_state.mid_g = g(base_state.mid_x, params);
base_state.mid_dcdz = dcdz(base_state.mid_x, params);
end

function out = block_matrix(segm_i, sigma, base_state, eqN)
% [Psi, X, Phi, Gamma]
%% set params
% def: C1(t) == sqrt(2t)/zeta2
C1 = base_state.params.C1;

R = base_state.params.R;
h = base_state.params.h;
h2 = h*h;

c = base_state.mid_c(segm_i);
g = base_state.mid_g(segm_i);
xBar = base_state.mid_xBar(segm_i);

%% set out
out = zeros(eqN,eqN);
out(1,1) = -g;
out(1,2) = -g*(1-c)./xBar;
out(1,3) = -R*g*(1-c);

out(2,1) = 1;
out(2,2) = (1-c)./xBar+xBar*(1+sigma);

out(3,4) = 1;

out(4,1) = h2;
out(4,3) = h2;
out(4,4) = -R*g/C1*(1-c)./xBar;
end

function out = Diag(segm_i, eqN, base_state)
% def: C1(t) == sqrt(2t)/zeta2
C1 = base_state.params.C1;
% def: C2(t) == 2t*dzeta2dt/zeta2
C2 = base_state.params.C2;

xBar = base_state.mid_xBar(segm_i);
xiMid = base_state.mid_xi(segm_i);

out = sparse(1:eqN, 1:eqN, ...
    [xBar*C1, xBar*C2*xiMid, 1, 1], eqN, eqN, eqN);
end

function out = BC_Keller(params, eqN)
% [Psi, X, Phi, Gamma]
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
h = params.h;

left = zeros(eqN,eqN);
right = zeros(eqN,eqN);
out.rhs = [0; 0; 0; -1];
%% Psi(0) = 0;
left(1,1) = 1;
%% Phi(0) = 0;
left(2,3) = 1;
%% Gamma(1)+h*Phi(1) = 0
right(3,3) = h;
right(3,4) = 1;
%% Gamma(1) = -1 --- set the scale
right(4,4) = 1;
%% set data
out.left = left;
out.right = right;
end

function out = BC_L(params, base_state, sigma, eqN, mesh)
% [Psi, X, Phi, Gamma]
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
[~,~, Psi_left, X_left, Phi_left, Gamma_left] = calc_full_inlet_solution_asymptotics(...
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
out.rhs = [Psi_left; X_left; Phi_left; Gamma_left];
end

function out = BC_R(params, base_state, sigma, eqN, mesh)
error('The method is not set')
% [Psi, X, Phi, Gamma]
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
% [Psi, X, Phi, Gamma]

% left*y(left) + right*y(right) = rhs
% [Psi] + a*dz0dt*g0*X = 0, Y = xi0*X
% - 1*Psi(left)   + 1*Psi(right) + a*dz0dt*g0*X(right) = 0
% - 1*X(left)     + 1*X(right) = 0
% - 1*Phi(left)   + 1*Phi(right) = 0
% - 1*Gamma(left) + 1*Gamma(right) = 0
a = params.a; % a == sqrt(2*t)
g0 = params.g0-params.g1;
dz2dt = params.dz2dt;
out.left = -eye(eqN,eqN);
out.right = eye(eqN,eqN);
out.right(1,2) = a*dz2dt*g0;
out.rhs = [0;0;0;0];
end

function sol = transform(sol, params, sigma, base_state)
% global xBarRight
t = reshape(sol.t, numel(sol.t), 1);% sol.t';
a0 = params.a0;
a = params.a;
dz0dt = params.dz0dt;
g1 = params.g1;
x_right = a0;

%% normalize
sol.factor = 1;
sol.y = sol.y/sol.factor;
%% input
% [Psi, X, Phi, Gamma]
Psi = sol.y(:,1);
X = sol.y(:,2);
Phi = sol.y(:,3);
Gamma = sol.y(:,4);
%% output
%[Psi, X, Phi, Gamma, Omega, Y, Psi+X]
Y = X.*t;
% Omega
g0 = params.g0;
dz2dt = params.dz2dt;
Omega = Psi + Y*a*g0*dz2dt;
%% flux Q, pressure P
Q = (Phi+Psi)./(1-t);
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
sol.y = [Psi, X, Phi, Gamma, Omega, Y, (Psi+X), Q, P, dP];
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
    -(g1*a*dz0dt + (1+sigma)*a0/a)*X(sol.id);
sol.sigma = sigma;
end

function dy = my_ode(x,y, sigma, params)
% [Psi, X, Phi, Gamma]
a = params.a; % == sqrt(2tau)
z = z_of_x(x, params);
C1 = params.C1;
C2 = params.C2;
z2 = params.z2;
h = params.h;
R = params.R;

gx = g(x, params);
Gx = G_of_x(x, params);
Ga = G_of_x(a, params);

f = zeros(4,4);
f(1,1) = gx/Gx;
f(1,2) = gx*a/(x*Ga);
f(1,3) = gx/Ga;

f(2,1) = -a/(Gx*z*C2);
f(2,2) = -(a*a/(x*Ga) + x*(1+sigma)/Gx)/(z*C2);

f(3,4) = -x/(z2*Gx);
h2 = h*h;

f(4,1) = -h2*R/(z2*G_div_X(x,params));
f(4,3) = -h2/(z2*G_div_X(x,params));
f(4,4) = R*gx/Ga;

dy = f*y;
end


