function u = IVBP_solver(u0, params)

t = params.t_grid; %linspace(1e-5,T, 1 + round(omega*(Nt-1)*T/period));

z = params.z;
m = params.m;
%% base state solution
x = x_of_z(z, t, params);
dcdz = zeros(size(x));
for i = 1:numel(t)
    dcdz(:,i) = dCdz(x(:,i), t(i), params);
end
params.x = x;
params.dcdz = dcdz;
%% initial-value calculation
u = zeros(2*m, numel(t));
u(:,1) = u0;
for i = 2:numel(t)
    params.iter = i;
%     i
    % current time step
    tau = t(i) - t(i-1);
    params.t = t(i);
    % matrix
    A = Matrix(tau, z, params);
    % rhs
    b = RHS(z, tau, u(:,i-1), params);
    % solve problem
    u(:,i) = A\b;
    % accumulate the rhs non-local term
    params.history = params.history + ...
        tau/2*(u(1:m,i-1).*dcdz(:, i-1) + u(1:m,i).*dcdz(:, i));    
end
end