function [h, out] = z_grid(h, L)
% nmbr of segments
n = ceil(L/h);
% z-grid mesh
out = linspace(0, L, n+1);
% z-grid uniform step
h = out(2) - out(1);
end