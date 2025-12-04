function sol = starter_R_zero_X_Psi(sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
% sol = transform(solver(problem, mesh), params, problem.base_state, sigma);
sol = transform(...
    solver_RK_backwards(problem.base_state, mesh, params, sigma), ...
    params, problem.base_state, sigma, mesh);
end

function problem = set_problem(mesh, sigma, params)
% number of equations
problem.eqN = 2;
problem.ids = 1:problem.eqN; % iterator for block rows/cols

problem.base_state = set_base_state(mesh, params);

problem.M = problem.eqN * mesh.N;
problem.block_matrix = @(xL, xR, i)block(xL, xR, i, sigma, params, problem.base_state, problem.eqN);
problem.diag = @(i)Diag(i, problem.eqN, problem.base_state);
problem.BC = @()BC(params, problem.base_state, sigma, problem.eqN, mesh);
problem.JC = @()JC(params, problem.base_state, problem.eqN);
end

function out = set_base_state(mesh, params)
a = params.a;

I = 1:(mesh.N-1);
nodes_xBar = [mesh.left.xBar, mesh.right.xBar];
nodes_xBar = (nodes_xBar(I)+nodes_xBar(I+1))/2;

out.z2 = z2(params);
out.z0 = z0(params);
out.dz2dt = dz2dt(params);
out.dz0dt = dz0dt(params);
out.C1 = a/out.z2;
out.C2 = out.C1*a*out.dz2dt;

out.x = nodes_xBar*a;
out.c = c_of_x(out.x, params);
out.g_of_x = g(out.x, params);
out.dcdz = dcdz(out.x, params);
out.dxBardt = dxBardt(out.x, params);
end

function out = block(xiL, xiR, i, sigma, params, base_state, eqN)
%[Psi, X]
% if(nargin == 3)
%     sigma = params.sigma;
% end

a = params.a;
%% set params
% def: C1(t) == sqrt(2t)/zeta2
C1 = base_state.C1;
% def: C2(t) == 2t*dzeta2dt/zeta2
C2 = base_state.C2;

c = base_state.c(i);
g_of_x = base_state.g_of_x(i);
xBar = base_state.x(i)/a;

y  =(xiL+xiR)/2;
%% set out
out = zeros(eqN,eqN);
out(1,1) = -g_of_x/C1/xBar;
out(1,2) = -(g_of_x/C1*(1-c)/(xBar^2));

out(2,1) = 1/(C2*xBar*y);
out(2,2) = ((1-c)/(xBar^2)+1+sigma)/y/C2;
end

function out = Diag(~, eqN, ~)
out = sparse(1:eqN, 1:eqN, ones(1,eqN), eqN, eqN, eqN);
end

function out = BC(params, base_state, sigma, eqN, mesh)
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
C2 = base_state.C2;
C1 = base_state.C1;
zeta2 = base_state.z2;
r = 1+(2+sigma)/C2;
a = params.a;
tau = (a*a)/2.0;
alpha = params.r;

xi_left = mesh.left.L;

[Psi_left, X_left] = ...
    calculate_X_Psi(xi_left, tau, alpha, C1, C2, zeta2, sigma);

%% BC at the left end
left = zeros(eqN,eqN);
left(1,1) = 1; % Psi(left) = Psi_left, left = delta -> 0
left(2,2) = 1; % X(left) = X_left, left = dalta -> 0
%% BC at the right end
right = zeros(eqN,eqN);
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
% X_left = Y_left/xi_left;
out.rhs = [Psi_left;X_left];
end

function out = JC(params, base_state, eqN)
%[Psi, Y]

% left*y(left) + right*y(right) = rhs
% [Psi] + a*dz0dt*g0*X = 0, Y = xi0*X
% - 1*Psi(left) + 1*Psi(right) + a*dz0dt*g0*X(right) = 0
% - 1*X(left) + 1*X(right) = 0
a = params.a; % a == sqrt(2*t)
%% JC at the left end
left = -eye(eqN,eqN);
%% rhs for BC eqns
% assign JC at the left end
out.left = left;
% JC at the right end
out.right = eye(eqN,eqN);
out.right(1,2) = a*base_state.dz0dt*(params.g0-params.g1);
out.rhs = [0;0];
end

function sol = transform(sol, params, base_state, sigma, mesh)
t = sol.t';
%% in
%[Psi, X]
sol.y = sol.y/sol.y(end,2);
Psi = sol.y(:,1);
X = sol.y(:,2);
Phi = 0*Psi;
Gamma = 0*Psi;
%% out
%[Phi, Psi, X, Gamma, Omega, Y, Psi+X]
% X
Y = X.*t;
% Omega
g0 = params.g0;
dz2dt = base_state.dz2dt;
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
sol.y = [Phi, Psi, X, Gamma, Omega, Y, 0*(Psi+X), 0*Q, 0*P, 0*dP];

%% Keller-box scheme behind the jump point
% C1 = base_state.C1;
% C2 = base_state.C2;
% sol.BC = Psi(end) + X(end)*g0/C1*(C2+(1+sigma)*(1-t(end)));
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
a = params.a; % == sqrt(2tau)
tau = (a*a)/2.0;
alpha = params.r;
C2 = base_state.C2;
C1 = base_state.C1;
zeta2 = base_state.z2;
xi_left = mesh.left.L;
[Psi_left, X_left] = ...
    calculate_X_Psi(xi_left, tau, alpha, C1, C2, zeta2, sigma);
sol.BC = Psi_left/X_left - Psi(1)/X(1);
sol.base_state = base_state;
sol.sigma = sigma;
end

function sol = solver_RK(base_state, mesh, params, sigma)
%% assymptotics at xi = 0
a = params.a; % == sqrt(2tau)
tau = (a*a)/2.0;
alpha = params.r;
C2 = base_state.C2;
C1 = base_state.C1;
zeta2 = base_state.z2;
xi_left = mesh.left.L;
[Psi_left, X_left] = ...
    calculate_X_Psi(xi_left, tau, alpha, C1, C2, zeta2, sigma);
%% init the RK solver
options = odeset('Abstol', 1e-10, 'RelTol', 1e-10);
y0 = [Psi_left; X_left];
x_mesh = mesh.left.xBar*a;
[t,y] = ode15s(@(x,y) my_ode(x,y,sigma,params,base_state), ...
    x_mesh, y0, options);
sol.t = z_of_x(t', params)/zeta2;
sol.y = y;
end

function sol = solver_RK_backwards(base_state, mesh, params, sigma)
%% assymptotics at xi = 0
a = params.a; % == sqrt(2tau)
zeta2 = base_state.z2;
%% init the RK solver
C2 = base_state.C2;
C1 = base_state.C1;
dz0dt = base_state.dz0dt;
g0 = params.g0 - params.g1;
g1 = params.g1;
xi0 = base_state.z0/base_state.z2;

options = odeset(...
    'Abstol', 1e-10...
    , 'RelTol', 1e-10 ...
    , 'NormControl', 'on' ...
    , 'NonNegative', 2 ...
    , 'MaxStep', 0.01);
X_right = 1;
Psi_right = (g0*a*dz0dt - (g0+g1)/C1*(C2+(1+sigma)*(1-xi0)));
y0 = [Psi_right; X_right];
x_mesh = mesh.left.xBar(end:-1:1)*a;
[t,y] = ode45(@(x,y) my_ode(x,y,sigma,params,base_state), ...
    x_mesh, y0, options);
sol.t = z_of_x(t', params)/zeta2;

sol.t = sol.t(end:-1:1);
sol.y = y(end:-1:1,:);
end

function dy = my_ode(x,y, sigma, params, base_state)

alpha = params.r;
a = params.a; % == sqrt(2tau)
g1 = 1-alpha;
D =  alpha+(1-alpha)*a;
Dx = alpha+(1-alpha)*x;
z = z_of_x(x, params);
C2 = base_state.C2;

f = zeros(2,2);
f(1,1) = g1/Dx;
f(1,2) = g1*a/(x*D);

f(2,1) = -a/(Dx*z*C2);
f(2,2) = -(a*a/(x*D) + x*(1+sigma)/Dx)/(z*C2);

dy = f*y;
end


