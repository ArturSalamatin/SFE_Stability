function func

clc
close all
clear all


params.t = 0.3;
params.a = sqrt(2*params.t);

params.R = 0.1;
params.alpha = 0.1;
sigma = 0.1;


opts = bvpset(...
    'RelTol',1e-5 ...
    , 'AbsTol',1e-5 ...
    , 'FJacobian',@(x,y) Jac_exp(x, y, sigma, params) ...
    , 'Stats','on');

A = 25;
K = 150;

xmesh = linspace(-A, A, 1001);
solinit = bvpinit(xmesh, @(t)left_expansion(t, K, sigma, params));

t = solinit.x;
ksi = exp(t)./(1+exp(t));
plot(ksi, solinit.y)

sol = bvp5c(...
    @(x,y) bvpfcn(x, y, params, sigma), ...
    @(ya, yb) bcfcn(ya, yb), ...
    solinit, opts);

plot(sol.x, sol.y(1,:), 'r-')

end

function dydx = bvpfcn(x, y, params, sigma)
dydx = Jac_exp(x, [], sigma, params)*y;
end

function res = bcfcn(ya, ~)
res = ya - [0;0;0;1];
end

function out = left_expansion(t, K, sigma, params)
gamma0 = 1;
R = params.R;

t = exp(t)./(1+exp(t));

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
