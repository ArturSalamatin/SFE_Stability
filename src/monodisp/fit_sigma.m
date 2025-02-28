function [sigma, sol] = fit_sigma(starter, guess, params)
sigma = guess;
[sigma, val] = ...
    fzero(@(s) func_to_min(starter, s, params), guess);
if(abs(val) > 1e-10)
    sigma = NaN;
%     error('Wrong value found!');
end
[val]
if(nargout == 2)
    sol = starter(sigma, params);
    plot_solution(sol.t,sol.y,params.pen)
end
end

function [out, sol] = func_to_min(starter, sigma, params)
sol = starter(sigma, params);
out = sol.y(end,5);
end

% function problem = set_problem(mesh, sigma, params)
% % number of equations
% problem.eqN = 4;
% problem.ids = 1:problem.eqN; % iterator for block rows/cols
% 
% problem.M = problem.eqN * mesh.N;
% problem.block_matrix = @(xL, xR)block_omega_Y(xL, xR, sigma, params);
% problem.BC = @()BC(params, problem.eqN);
% end
% 
% 
% function out = block_omega_Y(xiL, xiR, sigma, params)
% if(nargin == 3)
%     sigma = params.sigma;
% end
% 
% a = params.a;
% R = params.R;
% alpha2 = (params.alpha)^2;
% f2 = a*a*alpha2;
% 
% y  =(xiL+xiR)/2;
% 
% out = zeros(4,4);
% 
% out(1,4) = 1;
% 
% out(2,1) = -1;
% out(2,3) = (2+sigma)/y;
% 
% out(3,2) = 1/(1-y);
% out(3,3) = (3+sigma)/y;
% 
% out(4,1) = f2;
% out(4,2) = f2*R;
% out(4,3) = -f2*R;
% out(4,4) = -R;
% end
% 
% function out = BC_omega_Y(params, eqN)
% f = params.a*params.alpha;
% %% BC at the left end
% left = zeros(eqN,eqN);
% left(1,1) = 1; % Phi(0) = 0
% left(2,2) = 1; % Omega(0) = 0 /* = Psi(0)*/
% % left(3,3) = 1; % Y(0) = 0
% % left(4,4) = 1; % G(0)   = 1
% %% BC at the right end
% right = zeros(eqN,eqN);
% % right(2,2) = 1; % Omega(1) = 0
% right(3,1) = f; % f*Phi(1) + G(1) = 0
% right(3,4) = 1; 
% right(4,4) = 1; % G(1)   = 1
% %% rhs for BC eqns
% out.left = left;
% out.right = right;
% out.rhs = [0;0;0;1];
% end
% 
% function out = block_inf_domain(xiL, xiR, sigma, params)
% if(nargin == 3)
%     sigma = params.sigma;
% end
% 
% a = params.a;
% R = params.R;
% alpha2 = (params.alpha)^2;
% f2 = a*a*alpha2;
% 
% % -inf < y < +inf
% y  =(xiL+xiR)/2;
% % 0 < xi < 1
% xi = exp(y)/(1+exp(y));
% 
% out = zeros(4,4);
% 
% out(1,4) = xi*(1-xi);
% 
% out(2,1) = -xi*(1-xi);
% out(2,3) = xi*(1-xi)*(2+sigma);
% 
% out(3,2) = 1;
% out(3,3) = (1-xi)*(2+sigma);
% 
% out(4,1) = xi*(1-xi)*f2;
% out(4,2) = xi*(1-xi)*f2*R;
% out(4,3) = -xi*(1-xi)*f2*R*y;
% out(4,4) = -xi*(1-xi)*R;
% end
