function out = dGdx(x, params)
a0 = params.a0;
a1 = params.a1;
g0 = params.g0;
g1 = params.g1;
r = params.r;

out = ones(size(x));

mask = x <= a0;
out(mask) = g0;

mask = (x > a0) & (x < a1);
out(mask) = g1;

mask = x >= a1;
out(mask) = 0; % == (g0-g1)*a0+g1*(a1-a0);

out = out(:);
end