function out = z2(t, params)

a = sqrt(2*t);
a0 = params.a0;
g0 = params.g0;

out = min(a, a0)/g0 + z0(t, params);

end