function sol = solver_RK(problem, mesh, params)

%% assymptotics at xi = 0
bc = problem.BC_L();
y0 = bc.rhs;

% sigma = problem.sigma;
% g1 = params.g1;
% g0 = params.g0 - g1;
% xi0 = params.z0/params.z2;
% C1 = params.C1;
% C2 = params.C2;
% psi0 = -(g0+g1)/C1*(C2+(1+sigma)*(1-xi0));
% y0 = [psi0;1];

%% init the RK solver
options = odeset(...
    'Abstol', 1e-9...
    , 'RelTol', 1e-9 ...
    , 'NormControl', 'on' ...
    , 'MaxStep', 0.1);

[t,y] = ode15s(@(x,y) problem.RK(x,y), ...
    mesh.x, y0, options);

sol.x = t';
sol.t = z_of_x(t', params)/params.z2;
sol.y = y;
sol.id = numel(t);

jump = problem.JC();
y0 = -jump.right\(jump.left*sol.y(end,:)');
sol.y = sol.y/y0(2);
y0 = y0/y0(2);

if(mesh.x(end) < 0)
    error('Wrong mesh!');
end

[t2,y2] = ode45(@(x,y) problem.RK(x,y), ...
    [mesh.x(end)-1e-10, params.a0/5], y0, options);

sol.x = [sol.x, t'];
sol.t = [sol.t, z_of_x(t2', params)/params.z2];
sol.y = [sol.y;y2];
end