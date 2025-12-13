function sol = solver_RK(problem, mesh, params)
%% assymptotics at xi = 0
bc = problem.BC_L();

%% init the RK solver
options = odeset(...
    'Abstol', 1e-12...
    , 'RelTol', 1e-12 ...
    , 'NormControl', 'on' ...
    , 'NonNegative', 2 ...
    , 'MaxStep', 0.001);

[t,y] = ode45(@(x,y) problem.RK(x,y), ...
    mesh.x, bc.rhs, options);

jump = problem.JC();
y0 = -jump.right\(jump.left*y(end,:)');

if(mesh.x(end) < 0)
    error('Wrong mesh!');
end

[t2,y2] = ode15s(@(x,y) problem.RK(x,y), ...
    [mesh.x(end), 5e-2], y0, options);

sol.t = z_of_x([t;t2]', params)/params.z2;
sol.y = [y;y2];
sol.id = numel(t);
end