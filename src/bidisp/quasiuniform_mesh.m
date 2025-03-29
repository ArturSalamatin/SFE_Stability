function mesh = quasiuniform_mesh(L, R, N)
% mesh
mesh.L = L;
mesh.R = R;
mesh.N = N;
mesh.t = linspace(L, R, N);
mesh.I = 1:(N-1);
% steps
% mesh.steps = xi(I+1) - xi(I);
end