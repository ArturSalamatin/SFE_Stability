function out = dz2dt(params)

a = params.a;
a0 = params.a0;
r = params.r;
g0 = params.g0;

out = zeros(size(a));

% here z0 does not propagate yet
mask = a < a0;
out(mask) = 1./(a(mask)*g0);
% here z0 > 0, and front popagates
mask = ~mask;
out(mask) = 1./(r+(1-r)*a(mask));

end