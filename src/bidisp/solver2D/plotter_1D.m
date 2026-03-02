function [] = plotter_1D(states, mesh, params)

T = mesh.t(end);
H = mesh.z(end);

L = 1:(mesh.size)/mesh.Nr;
p = states(0*mesh.size+L,:);
c = states(1*mesh.size+L,:);
y = states(2*mesh.size+L,:);
G = states(3*mesh.size+L,:);
x = states(4*mesh.size+L,:);

figure(1005)
hold on
plot(mesh.z, x(:, end), 'k')
axis([0 H 0 1])
z = z_of_x(x(:,end), params, T);
plot(z, x(:, end), 'r')

figure(1002)
hold on
plot(mesh.z, c(:, end), 'k')
axis([0 H 0 1])
z = z_of_x(x(:,end), params, T);
C = c_of_x(x(:,end), params, T);
plot(z, C, 'r')

figure(1001)
hold on
plot(mesh.z, p(:, end), 'k')
axis([0 H 0 H])
% z = z_of_x(x(:,end), params, T);
% C = c_of_x(x(:,end), params, T);
% plot(z, C, 'r')

end

