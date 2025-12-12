function sol = solver_RK(problem, mesh, params)
%% assymptotics at xi = 0
bc = problem.BC_L();

%% init the RK solver
options = odeset('Abstol', 1e-10, 'RelTol', 1e-10);
[t,y] = ode15s(@(x,y) problem.RK(x,y), ...
    mesh.x, bc.rhs, options);

jump = problem.JC();
y0 = -jump.right\(jump.left*y(end,:)');

if(mesh.x(end) < 0)
    error('Wrong mesh!');
end

[t2,y2] = ode15s(@(x,y) problem.RK(x,y), ...
    [mesh.x(end), 1e-10], y0, options);

sol.t = z_of_x([t;t2]', params)/params.z2;
sol.y = [y;y2];
end