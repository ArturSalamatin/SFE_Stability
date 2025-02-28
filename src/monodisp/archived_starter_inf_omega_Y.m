function sol = starter_inf_omega_Y(sigma, params, mesh)
%% problem descriptor
problem = set_problem(mesh, sigma, params);
%% solve problem
sol = transform(solver(problem, mesh));
end

function problem = set_problem(mesh, sigma, params)
% number of equations
problem.eqN = 4;
problem.ids = 1:problem.eqN; % iterator for block rows/cols

problem.M = problem.eqN * mesh.N;
problem.block_matrix = @(xL, xR)block(xL, xR, sigma, params);
problem.BC = @()BC(params, problem.eqN);
end

function out = block(xiL, xiR, sigma, params)
if(nargin == 3)
    sigma = params.sigma;
end
% [Phi, Omega, Y, Gamma]

R = params.R;
f2 = (params.f)^2;

% -inf < y < +inf
y  =(xiL+xiR)/2;
% 0 < xi < 1
xi = exp(y)/(1+exp(y));

out = zeros(4,4);

out(1,4) = xi*(1-xi);

out(2,1) = -xi*(1-xi);
out(2,3) = (1-xi)*(2+sigma);

out(3,2) = xi;
out(3,3) = (1-xi)*(2+sigma);

out(4,1) = xi*(1-xi)*f2;
out(4,2) = xi*(1-xi)*f2*R;
out(4,3) = -xi*(1-xi)*f2*R;
out(4,4) = -xi*(1-xi)*R;
end

% function out = BC(params, eqN)
% f = params.f;
% %% BC at the left end
% left = zeros(eqN,eqN);
% left(1,1) = 1; % Phi(0) = 0
% left(2,2) = 1; % Omega(0) = 0 /* = Psi(0)*/
% % left(3,3) = 1; % Y(0) = 0
% left(4,4) = 1; % G(0)   = 1
% %% BC at the right end
% right = zeros(eqN,eqN);
% % right(2,2) = 1; % Omega(1) = 0
% right(3,1) = f; % f*Phi(1) + G(1) = 0
% right(3,4) = 1; 
% % right(4,4) = 1; % G(1)   = 1
% %% rhs for BC eqns
% out.left = left;
% out.right = right;
% out.rhs = [0;0;0;1];
% end

function sol = transform(sol)
t = sol.t';
xi = 1-1./(1+exp(t));
t = xi;
% in
% [Phi, Omega, Y, Gamma]
Phi = sol.y(:,1);
Omega = sol.y(:,2);
Y = sol.y(:,3);
Gamma = sol.y(:,4);
% out
%[Phi, Psi, X, Gamma, Omega, Y, Psi+X]
X = Y./t;
X(1) = 0;
Psi = Omega - Y;

sol.y = [Phi, Psi, X, Gamma, Omega, Y, Psi+X];
sol.t = t;
end