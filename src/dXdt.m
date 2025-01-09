function out = dXdt(X, params)
a = params.a;
out = 1/G_of_x(a)*G2X(X);
end