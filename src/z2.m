function out = z2(params)

a = params.a;
a0 = params.a0;
r = params.r;
g0 = params.g0;
g1 = params.g1;

if(a < a0)
    out = a/g0;
else
    out = a0/g0 + (a-a0+...
        r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*a)))/g1;
end
end