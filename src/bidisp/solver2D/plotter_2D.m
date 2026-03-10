function [] = plotter_2D(states, mesh)

Nr = mesh.Nr;
Nz = mesh.Nz;

L = 1:mesh.size;
p = states(0*mesh.size+L,:);
c = states(1*mesh.size+L,:);
y = states(2*mesh.size+L,:);
G = states(3*mesh.size+L,:);
x = states(4*mesh.size+L,:);

[r, z] = meshgrid(mesh.r, mesh.z);

Nt = mesh.Nt;
x_2D = reshape(x(:, Nt), Nz, Nr);
c_2D = reshape(c(:, Nt), Nz, Nr);

figure(2005)
hold on
contourf(r,z,x_2D)
figure(2002)
hold on
contourf(r,z,c_2D)

end

