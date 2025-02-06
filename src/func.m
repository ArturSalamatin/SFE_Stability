function func

clc
close all
clear all


params.t = 0.3;
params.a = sqrt(2*params.t);

params.R = 0.1;
params.alpha = 100;
sigma = 0.1;


opts = bvpset(...
    'RelTol',1e-5 ...
    , 'AbsTol',1e-5 ...
    , 'FJacobian',@(x,y) Jac(x, y, params, sigma) ...
    , 'Stats','on');

h = 1e-2;
x0 = h;
x1 = 1 - h;


xmesh = linspace(x0, x1, 10);
solinit = bvpinit(xmesh, [0; 0; 0; 1]);

sol = bvp4c(...
    @(x,y) bvpfcn(x, y, params, sigma), ...
    @(ya, yb) bcfcn(ya, yb, params, sigma, h), ...
    solinit, opts);

plot(sol.x, sol.y(1,:), 'r-')

end

function dydx = bvpfcn(x, y, params, sigma)
dydx = Jac(x, [], params, sigma)*y;
end

function res = bcfcn(ya, ~, params, sigma, h)
ic = IC(params, sigma, h);
res = ya - ic;
end

function yinit = mat4init(x) % initial guess function
yinit = [cos(4*x)
        -4*sin(4*x)];
end

function out = IC(params, sigma, h)
g0 = 1;
R = params.R;
alpha2 = (params.alpha)^2;
out = g0*[1*h - R/2*h*h
       -h*h/2
       h*h/(2*sigma)
       1-R*h+(R*R+alpha2)/2*h*h];
end

