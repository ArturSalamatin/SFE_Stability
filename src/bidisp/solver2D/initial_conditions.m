function IC = initial_conditions(mesh, params)
%INITIAL_CONDITIONS Summary of this function goes here
%   Detailed explanation goes here

IC.p = (params.H*ones(1,mesh.Nr) - mesh.z);
IC.p(end,:) = 0;
IC.c = ones(mesh.Nz, mesh.Nr);
IC.c(1,:) = 0;
IC.x = zeros(mesh.Nz, mesh.Nr)+1e-10;
IC.y = IC.x.*IC.x/2;
IC.G = params.g0*IC.x;
end
