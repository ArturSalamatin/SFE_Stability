function mesh = log_mesh(q, N, params)
% set mesh

xL = 1.0;
xR = 0.0;
% x0 == a0 <=> x0Bar == a0/sqrt(2t)
xMid = params.a0/params.a;

if(xMid > xL || xMid < xR)
    error('Wrong mesh stencils!');
end

NL = 1+ceil(N*(xL-xMid)/(xL-xR));
NR = 1+ceil(N*(xMid-xR)/(xL-xR));
N = NL+NR;

%% set left mesh segment
q = max(q, 1/q);
mesh.left = logMesh(...
    xL, xMid, q, NL, params);
%% set left mesh segment
mesh.right = logMesh(...
    xMid, xR, 1/q, NR, params);

mesh.N = mesh.right.N + mesh.left.N;
mesh.I = [mesh.left.I, mesh.left.N + mesh.right.I];

mesh.t = [mesh.left.t, mesh.right.t];
end

function mesh = logMesh(xL, xR, q, N, params)

h1 = (xR - xL)*(q-1)/(q^N-1);
i = 0:(N-1);
nodes = h1*q.^(i);
mesh.xBar = cumsum([xL, nodes]);
mesh.xBar([1,end]) = [xL, xR];

N = N+1;

% mesh
a = params.a;
mesh.xL = xL;
mesh.xR = xR;
mesh.N = N;
% mesh.xBar = linspace(xL, xR, N);
x = a*mesh.xBar;
mesh.t = z_of_x(x, params)/z2(params);
mesh.L = mesh.t(1);
mesh.R = mesh.t(end);
mesh.I = 1:(N-1);
end