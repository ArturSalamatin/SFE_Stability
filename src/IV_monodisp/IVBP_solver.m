function u = IVBP_solver(u0, params)

t = linspace(1e-5, 0.5, 2001);

z = params.z;
m = params.m;
%% base state solution
x = zeros(m, numel(t));
dcdz = zeros(size(x));
for i = 1:numel(t)
    x(:,i) = x_of_z(z, t(i));
    dcdz(:,i) = dCdz(x(:,i), t(i));
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