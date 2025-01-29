function out = dCdz(x, params)
a = params.a;
out = zeros(size(x));
mask = x > 0;
out(mask) = 1/a;
end