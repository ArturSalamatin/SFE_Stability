function out = z2(params)

a = params.a;
a0 = params.a0;
r = params.r;
g1 = params.g1;

if(a < a0)
    % here z0 does not propagate yet
    out = ...
        0;
else
    % here z0 > 0, and popagates
    out = (a-a0+...
        r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*a)))/g1;
end




end