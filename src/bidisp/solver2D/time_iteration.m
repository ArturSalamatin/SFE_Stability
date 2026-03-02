function state_s = time_iteration(...
    state_prev, mesh, params, dt)

state_s = state_prev;

p_id = 0*mesh.size + (1:mesh.size)';
c_id = 1*mesh.size + (1:mesh.size)';
y_id = 2*mesh.size + (1:mesh.size)';
G_id = 3*mesh.size + (1:mesh.size)';
x_id = 4*mesh.size + (1:mesh.size)';

s = 0;
while true
    d_s = single_iteration(...
        state_s, state_prev, mesh, params, dt);
    
    d_s(x_id) = max(d_s(x_id), -state_s(x_id));
    d_s(y_id) = max(d_s(y_id), -state_s(y_id));
    d_s(G_id) = max(d_s(G_id), -state_s(G_id));
    d_s(c_id) = max(d_s(c_id), -state_s(c_id));
    d_s(c_id) = min(d_s(c_id), 1-state_s(c_id));
    
    state_s = state_s + d_s;
    s = s+1;
    if(sum(abs(d_s)) < 1e-10)
        s
        break;
    end
end

end