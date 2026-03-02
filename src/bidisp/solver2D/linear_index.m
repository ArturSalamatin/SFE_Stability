function L = linear_index(J, I, mesh)
%LINEAR_INDEX Summary of this function goes here
%   Detailed explanation goes here

L = J + (I-1)*mesh.Nz;
L = L(:);
end

