function sol = solver_RK(problem, mesh, params)
%% assymptotics at xi = 0
bc = problem.BC_L();

%% init the RK solver
options = odeset('Abstol', 1e-10, 'RelTol', 1e-10);
[t,y] = ode15s(@(x,y) problem.RK(x,y), ...
    mesh.x, bc.rhs, options);
sol.t = z_of_x(t', params)/params.z2;
sol.y = y;
end