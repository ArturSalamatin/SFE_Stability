function IC = initial_conditions(mesh, params)
%INITIAL_CONDITIONS Summary of this function goes here
%   Detailed explanation goes here

global x_tol

IC.p = (params.H*ones(1,mesh.Nr) - mesh.z);
IC.p(end,:) = 0;
IC.c = ones(mesh.Nz, mesh.Nr);
IC.c(1,:) = 0;
IC.y = zeros(mesh.Nz, mesh.Nr)+x_tol*x_tol/2;
IC.x = sqrt(2*IC.y);
IC.G = params.g0*IC.x;
end
