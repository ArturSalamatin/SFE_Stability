function right_solver(params, sigma, fin, vals)

if(nargin == 0)
    clc
    % close all
    clear all
    
    sigma = 0.1;
    
    params.t = 0.1;
    params.a = sqrt(2*params.t);
    
    params.R = 0.2;
    params.alpha = 10;
    global A
    A = -35;
end

a = params.a;
a2 = a*a;
alpha = params.alpha;
alpha2 = alpha*alpha;
R = params.R;

if(nargin < 3)
    global A
    fin = 0.5;
    phi0 = 0.5;
    psi0 = 1;
    chi0 = -psi0;
    gamma0 = -a*alpha*phi0;
    
    IC = [phi0; psi0; chi0; gamma0];
    
    c1 = exp(A)/2;
    phi1 = -gamma0;
    chi1 = -phi0;
    psi1 = 2*phi0-(1+sigma)*chi0;
    gamma1 = R*gamma0 - a2*alpha2*(phi0+R*psi0);
    
    c2 = c1*c1;
    phi2 = -gamma1/2;
    chi2 = (gamma0 + (1+sigma)*(phi0-chi0))/4;
    psi2 = chi2-gamma0;
    gamma2 = (R*gamma1 - a2*alpha2*(phi1+R*psi1))/2;
    
    IC = IC + ...
        [phi1; psi1; chi1; gamma1]*c1 + ...
        [phi2; psi2; chi2; gamma2]*c2;
end
global fig_id A
fig_id = 700;

%% inf solver
Jac = @(z, u) Jac_plus(z, u, params, sigma);
[t, y] = ...
    solver(sigma, params, Jac, IC, [-A, 0]);
plot_solution(1-exp(-t)/2,y,'-r');

%% backward inf solver
[t, y] = ...
    solver(sigma, params, Jac, y(end,:), [0, -A]);
plot_solution(1-exp(-t)/2,y,'--k');
end

function dy = ode(t, y, sigma, params)

dy = J(t,y, sigma, params)*y;

end

function out = J(t, ~, sigma, params)

out = zeros(4,4);
out(1,4) = 1;

out(2,1) = -1;
out(2,3) = 1+sigma;

out(3,2) = 1;
out(3,3) = t + (1-t)*(2+sigma);

a = params.a;
alpha = params.alpha;
R = params.R;
out(4,1) = (a*alpha)^2;
out(4,2) = out(4,1)*R;
out(4,4) = -R;

end

function plot_solution(t,y,col)
global fig_id
for i = 1:4
    figure(fig_id+i)
    axis([0 1 -Inf Inf])
    plot(t,y(:,i),col, 'LineWidth', 1)
    hold on
end
end

function out = right_expansion(t, K, sigma, params)
gamma0 = 1;
R = params.R;


c0 = [0,0,0,gamma0];
c1 = [gamma0,0,0,-R*gamma0];

out = c0 + c1.*t;
c = c1;
T = t;
for i = 2:K
    c = coefs(c,i,params, sigma);
    T = T.*t;
    out = out + c.*T;
end

end

function out = coefs(c, k, params, sigma)
% [phi, psi, x, gamma]
R = params.R;
a = params.a;
alpha = params.alpha;

out = zeros(size(c));

out(1) = c(4)/(k+1);
out(2) = -(c(1) + (k-1-sigma)*c(3))/(k+1);
out(3) = out(2)/(k-1-sigma)+c(3);
out(4) = (((a*alpha)^2)*(c(1)+R*c(2)) - R*c(4))/(k+1);

end
