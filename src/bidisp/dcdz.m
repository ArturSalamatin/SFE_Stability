function out = dcdz(X, params)
a = params.a;
out = -g(X, params).*dXdz(X,params)/G_of_x(a, params);
end