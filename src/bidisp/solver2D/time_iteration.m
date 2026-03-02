function [state_s, s] = time_iteration(...
    state_prev, mesh, params, dt, t_prev)

global x_tol

Nz = mesh.Nz;

state_s = state_prev;

p_id = 0*mesh.size + (1:mesh.size)';
c_id = 1*mesh.size + (1:mesh.size)';
y_id = 2*mesh.size + (1:mesh.size)';
G_id = 3*mesh.size + (1:mesh.size)';
x_id = 4*mesh.size + (1:mesh.size)';

state_s(y_id) = state_s(y_id) + dt;
state_s(x_id) = sqrt(2*state_s(y_id));
ll = 3*mesh.size + (1:Nz:mesh.size)';
state_s(G_id) = G_of_x(state_s(x_id), params);
% state_s(G_id) = G_of_x(state_s(x_id), params);

s = 0;
while true
%     x = state_s(4*mesh.size + (1:mesh.size/2)');
%     G = state_s(3*mesh.size + (1:mesh.size/2)');
%     y = state_s(2*mesh.size + (1:mesh.size/2)');
%     c = state_s(1*mesh.size + (1:mesh.size/2)');
%     p = state_s(0*mesh.size + (1:mesh.size/2)');
    
    d_s = single_iteration(...
        state_s, state_prev, mesh, params, dt);
%     d_s(y_id) = ...
%         ... max(-state_s(y_id), ...
%         min(d_s(y_id), dt/5)...) ...
%         ;
    
%     d_x = d_s(4*mesh.size + (1:mesh.size/2)');
%     d_G = d_s(3*mesh.size + (1:mesh.size/2)');
%     d_y = d_s(2*mesh.size + (1:mesh.size/2)');
%     d_c = d_s(1*mesh.size + (1:mesh.size/2)');
%     d_p = d_s(0*mesh.size + (1:mesh.size/2)');
    
    d_s(x_id) = min(dt, max(d_s(x_id), -state_s(x_id)));
    d_s(y_id) = max(d_s(y_id), -state_s(y_id));
    d_s(G_id) = max(d_s(G_id), -state_s(G_id));
    d_s(c_id) = max(d_s(c_id), -state_s(c_id));
    d_s(c_id) = min(d_s(c_id), 1-state_s(c_id));
    
    
    state_s = state_s + d_s;
    
%     x = state_s(4*mesh.size + (1:mesh.size/2)');
%     G = state_s(3*mesh.size + (1:mesh.size/2)');
%     y = state_s(2*mesh.size + (1:mesh.size/2)');
%     c = state_s(1*mesh.size + (1:mesh.size/2)');
%     p = state_s(0*mesh.size + (1:mesh.size/2)');
    
    
    state_s(y_id) = max(x_tol*x_tol/2, state_s(y_id));
    state_s(y_id) = min(t_prev+dt, state_s(y_id));
    state_s(x_id) = sqrt(2*state_s(y_id));
    state_s(G_id) = G_of_x(state_s(x_id), params);
    state_s(G_id) = min(1, state_s(G_id));
    state_s(c_id) = max(0, state_s(c_id));
    state_s(c_id) = min(1, state_s(c_id));
    
%     x = state_s(4*mesh.size + (1:mesh.size/2)');
%     G = state_s(3*mesh.size + (1:mesh.size/2)');
%     y = state_s(2*mesh.size + (1:mesh.size/2)');
%     c = state_s(1*mesh.size + (1:mesh.size/2)');
%     p = state_s(0*mesh.size + (1:mesh.size/2)');
    
    s = s+1;
    if(sum(abs(d_s)) < 1e-10)
%         s
        break;
    end
end
    
end