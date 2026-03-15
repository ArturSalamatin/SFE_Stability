function out = c_of_x(x, params, t)
% calculate c(x),
% valid for any x >= 0
a = sqrt(2*t);

out = 1 - G_of_x(x, params)/G_of_x(a, params);
end