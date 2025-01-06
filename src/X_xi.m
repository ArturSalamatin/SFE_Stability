function [z, X] = X_xi(t, r)
global x_figid
xi = linspace(0, 1/(1-r)*log(G(sqrt(2*t), r)/r), 1001);

X = 1/(1-r)*(G(sqrt(2*t), r)*exp(-xi*(1-r)) - r);

z = 1/(1-r)*(sqrt(2*t) - X - r*xi);

% figure(x_figid)
% plot(z/z(end), X/sqrt(2*t), '-r')
