function [] = plotter_1D(states, mesh, params)

T = mesh.t(end);
H = mesh.z(end);

L = 1:(mesh.size)/mesh.Nr;
p = states(0*mesh.size+L,:);
c = states(1*mesh.size+L,:);
y = states(2*mesh.size+L,:);
G = states(3*mesh.size+L,:);
x = states(4*mesh.size+L,:);


L_wall = (1:mesh.Nz)+mesh.size - mesh.Nz;
p_wall = states(0*mesh.size+L_wall,:);
c_wall = states(1*mesh.size+L_wall,:);
y_wall = states(2*mesh.size+L_wall,:);
G_wall = states(3*mesh.size+L_wall,:);
x_wall = states(4*mesh.size+L_wall,:);

figure(1005)
hold on
plot(mesh.z, x(:, end), 'k')
plot(mesh.z, x_wall(:, end), 'b')
axis([0 H 0 1])
z = z_of_x(x(:,end), params, T);
plot(z, x(:, end), 'r')

figure(1002)
hold on
plot(mesh.z, c(:, end), 'k')
plot(mesh.z, c_wall(:, end), 'b')
axis([0 H 0 1])
z = z_of_x(x(:,end), params, T);
C = c_of_x(x(:,end), params, T);
plot(z, C, 'r')

figure(1001)
hold on
plot(mesh.z, p(:, end), 'k')
plot(mesh.z, p_wall(:, end), 'b')
axis([0 H 0 H])
% z = z_of_x(x(:,end), params, T);
% C = c_of_x(x(:,end), params, T);
% plot(z, C, 'r')

end

