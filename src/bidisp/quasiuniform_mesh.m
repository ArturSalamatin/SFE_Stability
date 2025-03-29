function mesh = quasiuniform_mesh(L, Mid, R, N)
% set mesh
% approximate quasi-uniform step
h = (R-L)/(N-1);
%% set left mesh segment
mesh.left = uniform_mesh(L, Mid, ceil((Mid - L)/h) +1);
%% set left mesh segment
mesh.right = uniform_mesh(Mid, R, ceil((R - Mid)/h) +1);

mesh.N = mesh.right.N + mesh.left.N;
end

function mesh = uniform_mesh(L, R, N)
% mesh
mesh.L = L;
mesh.R = R;
mesh.N = N;
mesh.t = linspace(L, R, N);
mesh.I = 1:(N-1);
end