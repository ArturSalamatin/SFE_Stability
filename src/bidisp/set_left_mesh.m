function mesh = set_left_mesh(N, params, delta)
% set interval endpoints
xBarL = 1.0;
xBarR = params.a0/params.a; % = a/sqrt(2t)

xBarL = min(1.0, xBarL);
if(xBarL < xBarR)
    % if delta is too big,
    % set xL somewhere between 1.0 and the jump point xMid
    xBarL = xBarR*0.1 + 1.0*0.9;    
end

if(xBarR > xBarL)
    error('Wrong mesh stencils!');
end
%% set left mesh segment
mesh = logMesh(...
    xBarL, xBarR, N ...ceil(abs(xMid - xL)/h) +1
    , params);

mesh.N = mesh.right.N + mesh.left.N;
mesh.I = [mesh.left.I, mesh.left.N + mesh.right.I];

mesh.t = [mesh.left.t, mesh.right.t];
end

function mesh = logMesh(xBarL, xBarR, N, params)
xBar = logMesh0(xBarL, xBarR, N); % just mesh nodes in xBar variable
N = numel(xBar);
a = params.a; % = sqrt(2t)
x = a*xBar; % xBar*sqrt(2t)

% assemble mesh struct
mesh.N = N;
mesh.xBar = xBar;
mesh.xBarL = xBar(1);
mesh.xBarR = xBar(end);
mesh.x = x;
mesh.xi = z_of_x(x, params)/z2(params);
mesh.xiL = mesh.xi(1);
mesh.xiR = mesh.xi(end);
mesh.segmIds = 1:(N-1); % ids of all segments
end

%% account for singular points at both ends,
% merges two meshes with log step increments from both ends
function xBar = logMesh0(xL, xR, N)
global q
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
% no point should enter the mesh twice
xBar = [xBar, xBar2(2:end)];
end

%% script for a uniform mesh
% h = abs(xBarR-xBarL)/(N-1);
% mesh.left = uniform_mesh(...
%     xL, xMid, N ...ceil(abs(xMid - xBarL)/h) +1
%     , params);
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
