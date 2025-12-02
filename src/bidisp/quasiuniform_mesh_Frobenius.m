function mesh = quasiuniform_mesh_Frobenius(delta, N, params)
% set mesh
q = 1.005;
xL = 1.0-delta;
xMid = params.a0/params.a;
xR = xMid*0.0;

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
% mesh.left = uniform_mesh(...
%     xL, xMid, N ...ceil(abs(xMid - xL)/h) +1
%     , params);
mesh.left = logMesh(...
    xL, xMid, q, N ...ceil(abs(xMid - xL)/h) +1
    , params);
%% set left mesh segment
mesh.right = uniform_mesh(...
    xMid, xR, 25, params);

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

function mesh = logMesh(xL, xR, q, N, params)
% h1 = (xR - xL)*(q-1)/(q^N-1);
% i = 0:(N-1);
% steps = h1*q.^(i);
mesh.xBar = logMesh0(xL, xR, q, N); % cumsum([xL, steps]);
mesh.xBar([1,end]) = [xL, xR];

N = numel(mesh.xBar);% N+1;

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

function xBar = logMesh0(xL, xR, q, N)
N = ceil(N/2);
xC = (xR+xL)/2;
h = (xC - xL)*(q-1)/(q^N-1);

i = 0:(N-1);
steps = h*q.^(i);
xBar = cumsum([xL, steps]);
xBar([1,end]) = [xL, xC];

i = (N-1):-1:0;
steps = h*q.^(i);
xBar2 = cumsum([xC, steps]);
xBar2([1,end]) = [xC, xR];
xBar = [xBar, xBar2(2:end)];
end