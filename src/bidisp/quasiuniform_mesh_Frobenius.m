function mesh = quasiuniform_mesh_Frobenius(delta, N, params)
% set mesh

xL = 1.0-delta;
xMid = params.a0/params.a;
xR = 0.0;

if(xL < xMid)
    % if delta is too big,
    % set xL somewhere between 1.0 and the jump point xMid
    xL = xMid*0.1 + 1.0*0.9;    
end

if(xMid > xL || xMid < xR)
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
mesh.I = [mesh.left.I, mesh.left.N + mesh.right.I];

mesh.t = [mesh.left.t, mesh.right.t];
end

function mesh = uniform_mesh(xL, xR, N, params)
% mesh
a = params.a;
mesh.xL = xL;
mesh.xR = xR;
mesh.N = N;
mesh.xBar = linspace(xL, xR, N);
x = a*mesh.xBar;
mesh.t = z_of_x(x, params)/z2(params);
mesh.L = mesh.t(1);
mesh.R = mesh.t(end);
mesh.I = 1:(N-1);
end