function out = C(t, X, params)
r = params.r;
a = params.a;
out = 1 - G_of_x(X, params)/G_of_x(a, params);
end