function out = c_of_x(x, params)
% calculate c(x),
% valid for any x >= 0
a = params.a;

out = 1 - G_of_x(x, params)/G_of_x(a, params);
end