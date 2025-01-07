function out = G_of_x(x, params)
a0 = params.a0;
g0 = params.g0;
g1 = params.g1;

out = ones(size(x));

mask = x <= a0;
out(mask) = g0*x(mask);

mask = (x > a0) & (x < 1);
out(mask) = g0*a0 + g1*(x(mask) - a0);

end