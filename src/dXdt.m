function out = dXdt(X, params)
a = params.a;
out = G2X(X)/G_of_x(a);
end