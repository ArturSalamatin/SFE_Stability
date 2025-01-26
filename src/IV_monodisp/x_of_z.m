function out = x_of_z(z, t)
[t, z] = meshgrid(t,z);
out = max(0,a(t)-z);
end
