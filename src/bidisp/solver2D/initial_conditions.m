function IC = initial_conditions(mesh, params)
%INITIAL_CONDITIONS Summary of this function goes here
%   Detailed explanation goes here

IC.p = zeros(mesh.Nz, mesh.Nr);
IC.c = zeros(mesh.Nz, mesh.Nr);
IC.x = zeros(mesh.Nz, mesh.Nr);
IC.y = IC.x.*IC.x/2;
IC.G = params.g0*IC.x;
end
