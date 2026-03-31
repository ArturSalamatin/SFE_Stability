function [state_s, s] = time_iteration(...
    state_prev, mesh, params, dt, t_prev)

global x_tol

state_s = state_prev;

c_id = 1*mesh.size + (1:mesh.size)';
y_id = 2*mesh.size + (1:mesh.size)';
G_id = 3*mesh.size + (1:mesh.size)';
x_id = 4*mesh.size + (1:mesh.size)';

state_s(y_id) = state_s(y_id) + dt*(1-state_s(c_id));
state_s(x_id) = sqrt(2*state_s(y_id));
state_s(G_id) = G_of_x(state_s(x_id), params);

eps = 1e-1;
% v_in = 1+eps*(2*rand(1, numel(mesh.I))-1);
v_in = ones(size(mesh.r))+eps*(2*rand(1, numel(mesh.I))-1);
%v_in(1:end/2) = 0.5;

p_out = 0*eps*(2*rand(1, numel(mesh.I))-1);

s = 0;
while true
    d_s = single_iteration(...
        state_s, state_prev, mesh, params, dt, v_in, p_out);
    d_s(x_id) = min(dt, max(d_s(x_id), -state_s(x_id)));
    d_s(y_id) = max(d_s(y_id), -state_s(y_id));
    d_s(G_id) = max(d_s(G_id), -state_s(G_id));
    d_s(c_id) = max(d_s(c_id), -state_s(c_id));
    d_s(c_id) = min(d_s(c_id), 1-state_s(c_id));
    
    state_s = state_s + d_s;
    
    state_s(y_id) = max(x_tol*x_tol/2, state_s(y_id));
    state_s(y_id) = min(t_prev+dt, state_s(y_id));
    state_s(x_id) = sqrt(2*state_s(y_id));
    state_s(G_id) = G_of_x(state_s(x_id), params);
    state_s(G_id) = min(1, state_s(G_id));
    state_s(c_id) = max(0, state_s(c_id));
    state_s(c_id) = min(1, state_s(c_id));
    
    s = s+1;
    if(sum(abs(d_s))/numel(d_s) < 1e-10)
        break;
    end
    if(s>1000)
        break;
    end
end
    
end