function out = G_of_x(x, params)
a0 = params.a0;
a1 = params.a1;
g0 = params.g0;
g1 = params.g1;
r = params.r;

out = ones(size(x));

mask = x <= a0;
out(mask) = g0*x(mask);

mask = (x > a0) & (x < a1);
out(mask) = r + g1*x(mask);

mask = x >= a1;
out(mask) = 1; % == (g0-g1)*a0+g1*(a1-a0);

end