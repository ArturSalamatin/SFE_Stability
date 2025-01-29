function out = x_of_z(z, params)
a = params.a;
z2 = params.z2;
out = a-z;

mask = z > z2;
out(mask) = 0;
end
