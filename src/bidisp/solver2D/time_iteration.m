function state_s = time_iteration(...
    state_prev, mesh, params, dt)

state_s = state_prev;

s = 0;
while true
    d_s = single_iteration(...
        state_s, state_prev, mesh, params, dt);
    state_s = state_s + d_s;
    s = s+1;
    if(sum(abs(d_s)) < 1e-10)
        s
        break;
    end
end

end