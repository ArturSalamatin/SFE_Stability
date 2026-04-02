% function keller
% clc
%
% global fig_id
% fig_id = 700;
%
% % clear_plots();
%
% %% set params
% params.t = 0.4;
% params.a = sqrt(2*params.t);
% params.R = 0.5;
% params.alpha = 1;
% sigma = -1.8774931801682;
% %% solve at physical domain
% sol = starter_orig(sigma, params);
% plot_solution(sol.t,sol.y,'-r')
% % A = full(A);
% sol.y(1,3)
% %% solve at infinite domain
% % sol = starter_inf(sigma, params);
% % plot_solution(sol.t,sol.y,'-r')
% % % A = full(A);
% % sol.y(1,3)
% end


function [sigma, sol] = fit_sigma(guess, params)
sigma = guess;
[sigma, val] = ...
    fzero(@(sigma) func_to_min(sigma, params), guess);
val
if(nargout == 2)
    sol = starter_orig(sigma, params);
    plot_solution(sol.t,sol.y,'-k')
end
end

function [out, sol] = func_to_min(sigma, params)
sol = starter_orig(sigma, params);
out = sol.y(end,2);
end


function sol = starter_orig(sigma, params)
%% set mesh
L = 0;
R = 1;
N = 5001;
mesh = set_mesh(L, R, N);
%% problem descriptor
problem = set_problem_orig_domain(mesh, sigma, params);
%% solve problem
sol = solver(problem, mesh, sigma, params);
end

function sol = starter_inf(sigma, params)
%% set mesh
A = 25;
L = -A;
R = A;
N = 5001;
mesh = set_mesh(L, R, N);
%% problem descriptor
problem = set_problem_inf_domain(mesh, sigma, params);
%% solve problem
sol = solver(problem, mesh, sigma, params);
sol.t = exp(sol.t)./(1+exp(sol.t));
end

function sol = solver(problem, mesh, sigma, params)
eqN = problem.eqN;
ids = problem.ids;
M = problem.M;
t = mesh.t;
I = mesh.I;
%% allocate memory
Diag = sparse(1:eqN, 1:eqN, ones(1,eqN), eqN, eqN, eqN);
A = spalloc(M,M,M*2);
%% fill in the matrix
block_pos = 0;
for i = I
    step = t(i+1) - t(i);
    block = block_orig_domain(t(i), t(i+1), sigma, params);
    % coef at y_i
    A(block_pos+ids, block_pos + ids) = ...
        block + Diag/step;
    % coef at y_(i+1)
    A(block_pos+ids, block_pos + ids + eqN) = ...
        block - Diag/step;
    % move to the next iteration
    block_pos = block_pos+eqN;
end
% BC
bc = problem.BC();
A(block_pos+ids, ids) = bc.left;
A(block_pos+ids, block_pos+ids) = bc.right;
%% rhs
b = zeros(M, 1);
b(block_pos+ids) = bc.rhs;
%% solution
sol.y = reshape(A\b, problem.eqN, mesh.N)';
sol.t = mesh.t;
end

function mesh = set_mesh(L, R, N)
% mesh
mesh.L = L;
mesh.R = R;
mesh.N = N;
mesh.t = linspace(L, R, N);
mesh.I = 1:(N-1);
% steps
% mesh.steps = xi(I+1) - xi(I);
end

function mesh = set_mesh_boundary_layer(L, R, N)
% mesh
mesh.L = L;
mesh.R = R;
t1 = linspace(L, L + 0.0002, N);
t2 = linspace(t1(end), R, N);
mesh.t = [t1, t2(2:end)];

mesh.N = numel(mesh.t);
mesh.I = 1:(mesh.N-1);
% steps
% mesh.steps = xi(I+1) - xi(I);
end

function problem = set_problem_orig_domain(mesh, sigma, params)
% number of equations
problem.eqN = 4;
problem.ids = 1:problem.eqN; % iterator for block rows/cols

problem.M = problem.eqN * mesh.N;
problem.block_matrix = @(xL, xR)block_orig_domain(xL, xR, sigma, params);
problem.BC = @()BC_orig_domain(params, problem.eqN);
end

function problem = set_problem_inf_domain(mesh, sigma, params)
% number of equations
problem.eqN = 4;
problem.ids = 1:problem.eqN; % iterator for block rows/cols

problem.M = problem.eqN * mesh.N;
problem.block_matrix = @(xL, xR)block_inf_domain(xL, xR, sigma, params);
problem.BC = @()BC_orig_domain(params, problem.eqN);
end

function plot_solution(t,y,col)
fig_id = 700;

for i = 1:4
    figure(fig_id+i)
    box on
    %     axis([0 1 -Inf Inf])
    hold on
    xlabel('{\xi}')
end

for i = 1:4
    figure(fig_id+i)
    box on
    %     axis([0 1 -Inf Inf])
    plot(t,-y(:,i)/y(end,4),col, 'LineWidth', 1)
    hold on
end

figure(fig_id+1)
ylabel('{\Phi}')

figure(fig_id+2)
ylabel('{\Omega}')

figure(fig_id+3)
ylabel(['{-\Psi}, ', 'X'])

figure(fig_id+4)
ylabel('{\Gamma}')

figure(fig_id+3)
plot(t,-(-y(:,2) + t'.*y(:,3))/y(end,4),['--' col(end)], 'LineWidth', 1)
end

function clear_plots()
global fig_id

for i = 1:5
    figure(fig_id+i)
    %     legend
    box on
    hold off
end

end

function out = block_orig_domain(xiL, xiR, sigma, params)
if(nargin == 3)
    sigma = params.sigma;
end

a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;
f2 = a*a*alpha2;

y  =(xiL+xiR)/2;

out = zeros(4,4);

out(1,4) = 1;

out(2,1) = -1;
out(2,3) = 2+sigma;

out(3,2) = 1/(y*(1-y));
out(3,3) = (2+sigma)/y;

out(4,1) = f2;
out(4,2) = f2*R;
out(4,3) = -f2*R*y;
out(4,4) = -R;
end

function out = BC_orig_domain(params, eqN)
f = params.a*params.alpha;
%% BC at the left end
left = zeros(eqN,eqN);
left(1,1) = 1; % Phi(0) = 0
left(2,2) = 1; % Omega(0) = 0 /* = Psi(0)*/
% left(3,3) = 1; % X(0) = 0
left(4,4) = 1; % G(0)   = 1
%% BC at the right end
right = zeros(eqN,eqN);
% right(2,2) = 1; % Omega(1) = 0
right(3,1) = f; % f*Phi(1) + G(1) = 0
right(3,4) = 1;
% right(4,4) = 1; % G(1)   = 1
%% rhs for BC eqns
out.left = left;
out.right = right;
out.rhs = [0;0;0;1];
end

function out = block_inf_domain(xiL, xiR, sigma, params)
if(nargin == 3)
    sigma = params.sigma;
end

a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;
f2 = a*a*alpha2;

% -inf < y < +inf
y  =(xiL+xiR)/2;
% 0 < xi < 1
xi = exp(y)/(1+exp(y));

out = zeros(4,4);

out(1,4) = xi*(1-xi);

out(2,1) = -xi*(1-xi);
out(2,3) = xi*(1-xi)*(2+sigma);

out(3,2) = 1;
out(3,3) = (1-xi)*(2+sigma);

out(4,1) = xi*(1-xi)*f2;
out(4,2) = xi*(1-xi)*f2*R;
out(4,3) = -xi*(1-xi)*f2*R*y;
out(4,4) = -xi*(1-xi)*R;
end
