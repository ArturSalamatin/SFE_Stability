function mesh = set_mesh(Nr, Nz, dt, T, params)
%SET_MESH Summary of this function goes here
%   Detailed explanation goes here

H = params.H;

mesh.Nr = Nr;
mesh.Nz = Nz;
Nt = ceil(T/dt)+3;
mesh.Nt = Nt;

mesh.I = 1:Nr;
mesh.I_l = 1:(Nr-1);
mesh.I_r = 2:Nr;
mesh.J = [1:Nz]';
mesh.J_l = [1:(Nz-1)]';
mesh.J_r = [2:Nz]';

mesh.r = linspace(0, 1, Nr);
mesh.z = linspace(0, H, Nz)';
mesh.t = linspace(0,T, Nt);

mesh.dr = mesh.r(2)-mesh.r(1);
mesh.dz = mesh.z(2)-mesh.z(1);
mesh.dt = mesh.t(2) - mesh.t(1);

r_mid = [mesh.r(1), (mesh.r(mesh.I_l)+mesh.r(mesh.I_r))/2, mesh.r(end)];
mesh.r_area = pi*(r_mid(2:end).^2 - r_mid(1:end-1).^2);
z_mid = [mesh.z(1); (mesh.z(mesh.J_l)+mesh.z(mesh.J_r))/2; mesh.z(end)];
dz = (z_mid(2:end) - z_mid(1:end-1));
mesh.V = mesh.r_area.*dz;
% mesh.V([1,end],:) = mesh.V([1,end],:)/2;
% mesh.V(:,[1,end]) = mesh.V(:,[1,end])/2;
mesh.z_area = 2*pi*dz.*r_mid(2:end-1);

mesh.size = Nr*Nz;

mesh.D_z = mesh.z_area/mesh.dr;
mesh.D_r = mesh.r_area/mesh.dz;
mesh.DBB_z = (params.B*params.B)*mesh.D_z;
end

