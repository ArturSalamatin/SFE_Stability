function sol = starter_R_zero(sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
sol = transform(solver(problem, mesh), params, problem.base_state);
end

function problem = set_problem(mesh, sigma, params)
% number of equations
problem.eqN = 2;
problem.ids = 1:problem.eqN; % iterator for block rows/cols

problem.base_state = set_base_state(mesh, params);

problem.M = problem.eqN * mesh.N;
problem.block_matrix = @(xL, xR, i)block(xL, xR, i, sigma, params, problem.base_state);
problem.diag = @(i)Diag(i, problem.eqN, problem.base_state);
problem.BC = @()BC(params, problem.eqN);
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
out.g_of_x = g(out.x, params);
out.dcdz = dcdz(out.x, params);
out.dxBardt = dxBardt(out.x, params);

out.B21 = a*out.dcdz;
out.B33 = (a^3)./out.x.*out.dxBardt;
end

function out = block(xiL, xiR, i, sigma, params, base_state, eqN)
%[Psi, Y]
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
out(1,1) = g_of_x/C1/xBar;
out(1,2) = (g_of_x/C1*(1-c)/(xBar^2))/y;

out(2,1) = 1/(C2*xBar);
out(2,2) = (C2+(1-c)/(xBar^2)+1+sigma)/y;
end

function out = Diag(~, eqN, ~)
out = sparse(1:eqN, 1:eqN, ones(1,eqN), eqN, eqN, eqN);
end

function out = BC(params, base_state, sigma, eqN, d_zeta)
% d_zeta -- small value close to zeta = 0,
% it is used to cut the singular point zeta = 0
C2 = base_state.C2;
C1 = base_state.C1;
zeta2 = base_state.zeta2;
r = (2+sigma)/C2;
psi = -(1-alpha)/((r+1)*C1*zeta2)*(d_zeta^(r+1));
x = (1 + (A^2*gamma1/C2 - (1-alpha)/((r+1)*C2*C1*zeta2))*d_zeta)*(d_zeta^r);
d_xi = d_zeta/zeta2;
y = x*d_xi;
%% BC at the left end
left = zeros(eqN,eqN);
left(1,1) = 1; % Y(0) = 0
left(2,2) = 1; % Psi(0) = 0
%% BC at the right end
right = zeros(eqN,eqN);
%% rhs for BC eqns
% left*y(0) + right*y(1) = rhs
out.left = left;
out.right = right;
out.rhs = [y;psi];
end

function out = JC(params, base_state, eqN)
%% JC at the left end
left = eye(eqN,eqN);
left(2, 3) = base_state.dz0dt*(params.g0-params.g1); % [Psi] + dz0dt*g0*X = 0
%% rhs for BC eqns
% assign JC at the left end
out.left = left;
% JC at the right end
out.right = -eye(eqN,eqN);
out.rhs = [0;0;0;0];
end

function sol = transform(sol, params, base_state)
t = sol.t';
%% in
%[Phi, Psi, Y, Gamma]
Phi = sol.y(:,1);
Psi = sol.y(:,2);
Y = sol.y(:,3);
Gamma = sol.y(:,4);
%% out
%[Phi, Psi, X, Gamma, Omega, Y, Psi+X]
% X
X = Y./t;
%X(1) = 0;
% Omega
g0 = params.g0;
dz2dt = base_state.dz2dt;
a = params.a;
Omega = Psi + Y*a*g0*dz2dt;
%% flux Q, pressure P
Q = -(Phi+Psi)./(1-t);
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
end