function mesh = set_mesh(Nr, Nz, params)
%SET_MESH Summary of this function goes here
%   Detailed explanation goes here

H = params.H;

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

end

