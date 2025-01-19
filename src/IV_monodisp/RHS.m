function out = RHS(z, tau, u, params)
% u --- solution at previous time step
% must be initialized with init_value(..)
omega = params.omega;
t = params.t;
R = params.R;
alpha2 = (params.alpha)^2;

iter = params.iter;
x = params.x(:, iter);
dcdz = params.dcdz(:, iter);

m = numel(z);

out = zeros(size(u));

I1 = 1:m;
I2 = I1+m;

out(I1) = ...
    - x.*(dcdz*tau/2.*u(I1) + ...
    params.history);
out(I2) = R*alpha2/tau*u(I2);
% BC for phi
out(I2(1)) = sin(2*pi*t*omega);
out(I2(end)) = 0;
% BC for psi
out(I1(1)) = 0;

end