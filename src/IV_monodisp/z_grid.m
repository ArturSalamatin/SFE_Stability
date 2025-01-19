function [h, out] = z_grid(h, L)
n = ceil(L/h);
out = linspace(0, L, n+1);
h = out(2) - out(1);
end