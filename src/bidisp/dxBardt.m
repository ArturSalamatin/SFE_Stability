function out = dxBardt(X, params)
a = params.a;
out = (dXdt(X, params) - X/(a*a))/a;
end