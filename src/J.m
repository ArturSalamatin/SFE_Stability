function out = J(sigma, params)
[z, y] = solver(sigma, params);

out = log10( ...
    ... abs(z2 - z_end) + ...
    sum(abs(y(end, 2))));
end




