function out = g(x, params)

g0 = params.g0;
g1 = params.g1;
a0 = params.a0;
a1 = params.a1;

if(x < 0)
    error('X must be positive!');
end

if(x < a0)
    out = g0;
elseif(x < a1)
    out = g1;
else
    out = 0;
end
end