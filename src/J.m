function out = J(sigma, t, r, R, alpha, xi_max)

[~,y] = solver(sigma, t, r, R, alpha, xi_max);

out = log10(abs(y(end, 1)));

end

function [xi,y] = solver(sigma, t, r, R, alpha, xi_max)

options = odeset(...
      'RelTol', 1e-6 ...
    , 'AbsTol', 1e-6 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
    ..., 'InitialStep', 1e-8 ...
    ..., 'MaxStep', 1e-5 ...
    , 'Jacobian',@(xi, y) Jac(xi, y, t, r, R, alpha, sigma, xi_max) ...
    ..., 'Stats','on' ...
    ...,'OutputFcn', @odeplot ...
    );


[xi,y] = ode45(...
    @(xi, y) my_ode(xi, y, t, r, R, alpha, sigma, xi_max), ...
    [0, 1], [0;0;-1], options);

figure(701)
hold on
plot(xi, y(:,1))

figure(702)
hold on
plot(xi, y(:,2))

figure(703)
hold on
plot(xi, y(:,3))

[z, X] = X_xi(t, r);
c = C(t, X, r);
figure(200)
hold on
f = 1-c + sigma*X.^2;
plot(z/z(end), f)

% global X0 C0
% [X0, C0, xi(end), xi_max, sqrt(2*t)]
% Jac(xi(end), y(end,:), t, r, R, alpha, sigma)


end



