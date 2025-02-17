function f = left_solver(params, sigma)

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
    A = 35;
    
end
global A

fin = 0.9;

options = odeset(...
    'RelTol', 1e-14 ...
    , 'AbsTol', 1e-14 ...
    , 'NormControl', 'yes' ...
    , 'MaxStep', 1e-3 ...
    , 'InitialStep', 1e-5 ...
    , 'Jacobian', @(t,y) J(t,y,sigma, params) ...
    , 'Mass', @(t) mass(t) ...
    , 'MassSingular', 'yes' ...
    , 'MStateDependence', 'none' ...
    , 'InitialSlope', [1 0 0 -params.R] ...
    , 'BDF', 'on' ...
    );
%% DAE solver
[t,Y] = ode23t(...
    @(t,y) ode(t,y,sigma,params), ...
    [0,fin], ...
    [0,0,0,1], ...
    options);
x = fin;
u = Y(end,:);
plot_solution(t,Y,'-r', 'mass & ode23t 1e-14');
% [t,Y] = ode15s(...
%     @(t,y) ode(t,y,sigma,params), ...
%     [0,fin], ...
%     [0,0,0,1], ...
%     options);
% x = fin;
% u = Y(end,:);
% plot_solution(t,Y,'-k');
%% series expansion
K = 405;
t = linspace(0,0.9,301)';
out = left_expansion(t, K, sigma, params);
out(end, [1,4])
plot_solution(t, out, '--k', 'series expansion');
% K = 125;
% t = linspace(0,1,301)';
% out = left_expansion(t, K, sigma, params);
% plot_solution(t, out, '--r');
%% backward inf-solver
% Jac = @(z, u) Jac_minus(z, u, params, sigma);
% [t, y] = ...
%     solver(sigma, params, Jac, Y(end,:), [log(2*fin), -abs(A)]);
% plot_solution(exp(t)/2,y,'-b');
%% backward DAE solver
% [t,y] = ode23t(...
%     @(t,y) ode(t,y,sigma,params), ...
%     [fin,0.1], ...
%     Y(end,:), ...
%     options);
% plot_solution(t,y,'-r');
%% exp-vars
options = odeset(...
    'RelTol', 1e-12 ...
    , 'AbsTol', 1e-12 ...
    , 'NormControl', 'yes' ...
    ..., 'MaxStep', 1e-3 ...
    ..., 'InitialStep', 1e-5 ...
    , 'Jacobian', @(t,y) Jac_exp(t,y,sigma, params) ...
    ..., 'BDF', 'on' ...
    );

t = abs(35);
ksi = exp(-t)/(1+exp(-t));
out = left_expansion(ksi, 1, sigma, params);

[t,Y] = ode23t(...
    @(t,y) ode_exp(t,y,sigma,params), ...
    [-t,t], ...
    out, ...
    options);
x = t(end);
u = Y(end,:);
f = Y(end,1)/Y(end,4);
plot_solution(exp(t)./(1+exp(t)),Y,'-b','ode23t var_transform');
%% backward exp-vars
% [t,Y] = ode15s(...
%     @(t,y) ode_exp(t,y,sigma,params), ...
%     [x,-5], ...
%     u, ...
%     options);
% plot_solution(exp(t)./(1+exp(t)),Y,'--m');
end

function dy = ode(t, y, sigma, params)
dy = J(t,[], sigma, params)*y;
end

function dy = ode_exp(t, y, sigma, params)
dy = Jac_exp(t,[], sigma, params)*y;
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

function out = mass(t)

out = eye(4,4);
out(2,3) = t;
out(3,3) = t*(1-t);

end

function plot_solution(t,y,col,leg)
global fig_id
if(nargin == 3)
    leg = [];
end
for i = [1:4]
    figure(fig_id+i)
    axis([0 1 -Inf Inf])
    plot(t,y(:,i),col, 'LineWidth', 1, 'DisplayName', leg)
    hold on
end

    figure(fig_id+5)
    axis([0 1 -Inf Inf])
    plot(t,y(:,2)+y(:,3),col, 'LineWidth', 1, 'DisplayName', leg)
    hold on
end

function out = left_expansion(t, K, sigma, params)
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



