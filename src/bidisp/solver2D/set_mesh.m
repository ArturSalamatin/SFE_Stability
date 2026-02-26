function mesh = set_mesh(Nr, Nz, params)
%SET_MESH Summary of this function goes here
%   Detailed explanation goes here

H = params.H;

mesh.Nr = Nr;
mesh.Nz = Nz;

mesh.I = 1:Nr;
mesh.I_l = 1:(Nr-1);
mesh.I_r = 2:Nr;
mesh.J = [1:Nz]';
mesh.J_l = [1:(Nz-1)]';
mesh.J_r = [2:Nz]';

mesh.r = linspace(0, 1, Nr);
mesh.z = linspace(0, H, Nz)';

mesh.dr = mesh.r(2)-mesh.r(1);
mesh.dz = mesh.z(2)-mesh.z(1);

r_mid = [mesh.r(1), (mesh.r(mesh.I_l)+mesh.r(mesh.I_r))/2, mesh.r(end)];
mesh.r_area = pi*(r_mid(2:end).^2 - r_mid(1:end-1).^2);
z_mid = [mesh.z(1); (mesh.z(mesh.J_l)+mesh.z(mesh.J_r))/2; mesh.z(end)];
mesh.z_area = z_mid(2:end) - z_mid(1:end-1);
mesh.V = mesh.r_area.*mesh.z_area;

end

