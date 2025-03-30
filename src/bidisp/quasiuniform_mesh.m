function mesh = quasiuniform_mesh(xL, xMid, xR, N, params)
% set mesh
if(xL ~= 1.0 || xR ~= 0.0 || xMid > xL || xMid < xR)
    error('Wrong mesh stencils!');
end
% approximate quasi-uniform step in xBar(!)
h = abs(xR-xL)/(N-1);
%% set left mesh segment
mesh.left = uniform_mesh(...
    xL, xMid, ceil(abs(xMid - xL)/h) +1, params);
%% set left mesh segment
mesh.right = uniform_mesh(...
    xMid, xR, ceil(abs(xR - xMid)/h) +1, params);

mesh.N = mesh.right.N + mesh.left.N;
mesh.I = [mesh.left.I, mesh.left.N+1 + mesh.right.I];
end

function mesh = uniform_mesh(xL, xR, N, params)
% mesh
a = params.a;
mesh.xL = xL;
mesh.xR = xR;
mesh.N = N;
mesh.xBar = linspace(xL, xR, N);
mesh.t = z_of_x(a*mesh.xBar, params);
mesh.L = mesh.t(1);
mesh.R = mesh.t(end);
mesh.I = 1:(N-1);
end