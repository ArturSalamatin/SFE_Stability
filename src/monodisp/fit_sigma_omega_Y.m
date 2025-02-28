function [sigma, sol] = fit_sigma_omega_Y(guess, params)
sigma = guess;
[sigma, val] = ...
    fzero(@(s) func_to_min_omega_Y(s, params), guess);
[sigma, val]
if(nargout == 2)
    sol = starter_omega_Y(sigma, params);
    plot_solution(sol.t,sol.y,'-k')
end
end

function [out, sol] = func_to_min_omega_Y(sigma, params)
sol = starter_omega_Y(sigma, params);
out = sol.y(end,2);
end

function sol = starter_omega_Y(sigma, params)
%% set mesh
L = 0;
R = 1;
N = 5001;
mesh = set_mesh(L, R, N);
%% problem descriptor
problem = set_problem_omega_Y(mesh, sigma, params);
%% solve problem
sol = solver(problem, mesh);
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
