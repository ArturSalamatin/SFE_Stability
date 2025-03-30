function sol = starter_Y(sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
sol = transform(solver(problem, mesh), params, problem.base_state);
end

function problem = set_problem(mesh, sigma, params)
% number of equations
problem.eqN = 4;
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
out.C2 = a*a*out.dz2dt/out.z2;
out.C1 = a/out.z2;

out.x = nodes_xBar*a;
out.g_of_x = g(out.x, params);
out.dcdz = dcdz(out.x, params);
out.dxBardt = dxBardt(out.x, params);

out.B21 = a*out.dcdz;
out.B33 = (a^3)./out.x.*out.dxBardt;
end

function out = block(xiL, xiR, i, sigma, params, base_state)
%[Phi, Psi, Y, Gamma]
% if(nargin == 3)
%     sigma = params.sigma;
% end

a = params.a;
R = params.R;
h2 = (params.h)^2;

%% set params
% C2(t) = 2t*dzeta2dt/zeta2
C2 = base_state.C2;
% sqrt(2t)*dc/dzeta
B21 = base_state.B21(i);
% 2t/xBar*dxBardt
B33 = base_state.B33(i);

g_of_x = base_state.g_of_x(i);
dcdz = base_state.dcdz(i);
xBar = base_state.x(i)/a;


y  =(xiL+xiR)/2;
%% set out
out = zeros(4,4);

out(1,4) = 1;

out(2,1) = -B21;
out(2,3) = (1+sigma+C2)*g_of_x/y;

out(3,2) = 1/xBar/C2;
out(3,3) = (3+sigma+B33)/y/C2;

out(4,1) = h2;
out(4,2) = h2*R;
out(4,4) = -R*dcdz;
end

function out = Diag(i, eqN, base_state)
out = sparse([1:eqN,2], [1:eqN,3], ones(1,eqN+1), eqN, eqN, eqN+1);
out(2, 2) = base_state.C1;
out(2, 3) = base_state.C2*base_state.g_of_x(i);
end
function out = BC(params, eqN)
h = params.h;
%% BC at the left end
left = zeros(eqN,eqN);
left(1,1) = 1; % Phi(0) = 0
left(2,2) = 1; % Omega(0) = 0 /* = Psi(0)*/
% left(3,3) = 1; % Y(0) = 0
% left(4,4) = 1; % G(0)   = 1
%% BC at the right end
right = zeros(eqN,eqN);
% right(2,2) = 1; % Omega(1) = 0
right(3,1) = h; % f*Phi(1) + G(1) = 0
right(3,4) = 1; 
right(4,4) = 1; % G(1)   = 1
%% rhs for BC eqns
out.left = left;
out.right = right;
out.rhs = [0;0;0;-1];
end

function out = JC(params, base_state, eqN)
%% BC at the left end
left = eye(eqN,eqN);
left(2, 3) = base_state.dz0dt*(params.g0-params.g1); % [Psi] + dz0dt*g0*X = 0
%% rhs for BC eqns
% BC at the left end
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
X(1) = 0;
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