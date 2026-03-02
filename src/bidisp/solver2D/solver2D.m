function states = solver2D(mesh, params)

%% set initial conditions
IC = initial_conditions(mesh, params);
state = make_state(IC);

states = zeros(numel(state), numel(mesh.t));
states(:,1) = state;

for t_id = 2:numel(mesh.t)
    dt = mesh.dt;
    t_id*dt
    state = time_iteration(...
        state, mesh, params, dt, (t_id-1)*dt);
    
    states(:,t_id) = state;
end