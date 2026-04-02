function out = z2(params)

a = params.a;
a0 = params.a0;
g0 = params.g0;

out = min(a, a0)/g0 + z0(params);

end