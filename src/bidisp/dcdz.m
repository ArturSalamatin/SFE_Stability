function out = dcdz(X, params)
out = -g(X, params).*dXdz(X,params)/G_of_x(a, params);
end