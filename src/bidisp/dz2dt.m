function out = dz2dt(params)

a = params.a;
a0 = params.a0;
g0 = (params.g0 - params.g1);

if(a < a0)
    % here z0 does not propagate yet
    out = 1/(a*g0);
else
    % here z0 > 0, and front popagates
    out = 1/(r+(1-r)*a);
end
end