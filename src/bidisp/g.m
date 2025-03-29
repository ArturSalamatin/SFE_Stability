function out = g(x, params)

g0 = params.g0;
g1 = params.g1;
a0 = params.a0;
a1 = params.a1;

out = zeros(size(x));

for i = 1:numel(out)
if(x(i) < 0)
    error('X must be positive!');
end



if(x(i) < a0)
    out(i) = g0
elseif(x(i) == a0)
    out(i) = g0-g1;
elseif(x(i) < a1)
    out(i) = g1;
elseif(x(i) == a1)
    out(i) = g1 - 0;
else
    out(i) = 0;
end
end