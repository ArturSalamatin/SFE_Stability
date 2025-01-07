function out = dXdt(X, params)
a = params.a;
out = 1/G_of_x(a, params)*G2X(X, params);
end