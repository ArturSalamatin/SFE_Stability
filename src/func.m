function func

close all
clear all

x0 = 1/(30*pi);
options = odeset(...
    'RelTol', 1e-8 ...
    , 'AbsTol', 1e-8 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
     , 'InitialStep', 1e-8 ...
    , 'MaxStep', 1e-3 ...
    ... , 'Jacobian', Jac ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    ..., 'Events', @(z,y) events(z,y) ...
    );

[x, y] = ...
    ode45(@(x,y) f(x,y), ...
    [x0, 1], ...
    [sin(1/x0), -1/(x0*x0)*cos(1/x0)], options);
figure(1)
hold on
plot(x,y(:,1))


opts = bvpset(...
    'FJacobian',@jac,'RelTol',1e-5,...
    'AbsTol',1e-5,'Stats','on');
xmesh = linspace(x0, 1, 10);
solinit = bvpinit(xmesh, [1; 1]);

sol4c = bvp5c(@bvpfcn, @bcfcn, solinit, opts);

plot(sol4c.x,sol4c.y(1,:),'r*')

set(gca, 'xScale', 'log')

end



function dy = f(x, y)

dy = [y(2)
       -2*y(2)/x - y(1)/x^4];

end

function dfdy = jac(x,~)
dfdy = [0      1
       -1/x^4 -2/x];
end

function dydx = bvpfcn(x,y)
dydx = [y(2)
       -2*y(2)/x - y(1)/x^4];
end

function res = bcfcn(ya,yb)
res = [ya(1)
       yb(1)-sin(1)];
end

