function out = h0(z, t)

x = x_of_z(z,t);
[t, z] = meshgrid(t,z);
out = t - (x.^2)/2;


end