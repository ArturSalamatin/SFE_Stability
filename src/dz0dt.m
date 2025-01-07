function out = dz0dt(params)

a = params.a;
a0 = params.a0;
r = params.r;

if(a < a0)
    % here z0 does not propagate yet
    out = ...
        0;
else
    % here z0 > 0, and popagates
    out = 1/(r+(1-r)*a);
end



end