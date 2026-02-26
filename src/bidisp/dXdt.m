function out = dXdt(X, params)
a = params.a;
out = 1/G_of_x(a, params)*G_div_X(X, params);
end