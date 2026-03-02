function d_s = single_iteration(...
    state_s, state_prev, mesh, params, dt)
%SINGLE_ITERATION Summary of this function goes here
%   Detailed explanation goes here

% ratio of two spatial scales
B = params.B;
% nmbr of unknowns/equations
eq_nmbr = 5;
Nz = mesh.Nz;

p_id = 0*mesh.size + (1:mesh.size)';
c_id = 1*mesh.size + (1:mesh.size)';
y_id = 2*mesh.size + (1:mesh.size)';
G_id = 3*mesh.size + (1:mesh.size)';
x_id = 4*mesh.size + (1:mesh.size)';

I = mesh.I;
I_l = mesh.I_l;
I_r = mesh.I_r;
J = mesh.J;
J_l = mesh.J_l;
J_r = mesh.J_r;
rows = mesh.Nz;
cols = mesh.Nr;

p = reshape(state_s(p_id), rows, cols);
c = reshape(state_s(c_id), rows, cols);
y = reshape(state_s(y_id), rows, cols);
y_prev = reshape(state_prev(y_id), rows, cols);
G = reshape(state_s(G_id), rows, cols);
G_prev = reshape(state_prev(G_id), rows, cols);
x = reshape(state_s(x_id), rows, cols);
% p = state_s(p_id);
% c = state_s(c_id);
% y = state_s(y_id);
% y_prev = state_prev(y_id);
% G = state_s(G_id);
% G_prev = state_prev(G_id);
% x = state_s(x_id);

c_mid_z = (c(J_r,I)+c(J_l,I))/2;
c_mid_r = (c(J,I_r)+c(J,I_l))/2;
k_mid_z = K(c_mid_z, params);
k_mid_r = K(c_mid_r, params);
dk_mid_z = DK(c_mid_z, params);
dk_mid_r = DK(c_mid_r, params);

S_z = k_mid_r.*mesh.D_z;
SBB_z = B*B*S_z;
S_r = k_mid_z.*mesh.D_r;

q_z = -S_r.*(p(J_r,I)-p(J_l,I));
qB_r = -B*S_z.*(p(J,I_r)-p(J,I_l));
q_z_pos = (q_z + abs(q_z))/2;
q_z_neg = (q_z - abs(q_z))/2;
qB_r_pos = (qB_r + abs(qB_r))/2;
qB_r_neg = (qB_r - abs(qB_r))/2;

PBB_z = mesh.DBB_z.*(p(J,I_r)-p(J,I_l)).*dk_mid_r/2;
P_r = mesh.D_r.*(p(J_r,I)-p(J_l,I)).*dk_mid_z/2;

rhs = zeros(eq_nmbr*mesh.size, 1);
%% set rhs eqn: div(v*c) = dG/dt
shift = mesh.size;

L = linear_index(J, I, mesh);
rhs(L + shift) = mesh.V(:).*(G(:) - G_prev(:))/dt;

% j-1, i
L = linear_index(J_r, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - c(L-1).*q_z_pos(:);
% j, i
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*q_z_neg(:);
%j+1, i
L = linear_index(J_l, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) + c(L+1).*q_z_neg(:);
%j,i
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*q_z_pos(:);

%j, i
L = linear_index(J, I_r, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*qB_r_pos(:);
%j, i-1
rhs(L + shift) = ...
    rhs(L + shift) - c(L-Nz).*qB_r_neg(:);
%j, i
L = linear_index(J, I_l, mesh);
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*qB_r_neg(:);
%j, i+1
rhs(L + shift) = ...
    rhs(L + shift) + c(L+Nz).*qB_r_pos(:);
%% set rhs eqn: div(v) = 0
shift = 0;
% j-1, i
L = linear_index(J_r, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L-1).*S_r(:);
% j, i
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*S_r(:);
%j, i
L = linear_index(J_l, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L+1).*S_r(:);
%j+1,i
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*S_r(:);

%j, i
L = linear_index(J, I_r, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*SBB_z(:);
%j, i-1
rhs(L + shift) = ...
    rhs(L + shift) + p(L-Nz).*SBB_z(:);
%j, i
L = linear_index(J, I_l, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*SBB_z(:);
%j, i+1
rhs(L + shift) = ...
    rhs(L + shift) + p(L+Nz).*SBB_z(:);
%% set rhs eqn: dy/dt - (1-c) = 0
shift = 2*mesh.size;
L = linear_index(J, I, mesh);
rhs(L + shift) = -((y(:)-y_prev(:))/dt - (1-c(:)));
%% set rhs eqn: G(x) - ... = 0
shift = 3*mesh.size;
L = linear_index(J, I, mesh);
rhs(L + shift) = -(G - G_of_x(x, params));
%% set rhs eqn: y - x*x/2 = 0
shift = 4*mesh.size;
L = linear_index(J, I, mesh);
rhs(L + shift) = -(y(:) - x(:).*x(:)/2);
%% BC in rhs
% c_in = 0 in eq2
shift = mesh.size;
L = linear_index(J(1), I, mesh);
rhs(L + shift) = 0;
% v*c_out in eq2
shift = mesh.size;
L = linear_index(J(end), I, mesh);
rhs(L + shift) = rhs(L + shift) - c(L).*q_z_pos(J(end-1), I)';
% v_in = 1 in eq1
shift = 0;
L = linear_index(J(1), I, mesh);
rhs(L + shift) = rhs(L + shift) - mesh.D_r'*mesh.dz;
% p_out = 0
shift = 0;
L = linear_index(J(end), I, mesh);
rhs(L + shift) = 0;
%% problem matrix
m = eq_nmbr*mesh.size;
nnz = 7*mesh.size + 3*5*mesh.size;

L = linear_index(J, I, mesh)';

L1 = linear_index(J_l, I, mesh)';
L1_p = L1+1;
L2 = linear_index(J_r, I, mesh)';
L2_n = L2-1;
L3 = linear_index(J, I_l, mesh)';
L3_p = L3+Nz;
L4 = linear_index(J, I_r, mesh)';
L4_n = L4-Nz;

%rows
i = [... p-rows in eq1
    L1,L2,L3,L4, ...
    L1,L2,L3,L4, ...
    ... c-rows in eq1
    L1,L2,L3,L4, ...
    L1,L2,L3,L4, ...
    ... c-rows in eq2
    [L1,L2,L3,L4, ...
    L1,L2,L3,L4] + mesh.size, ...
    ... G-rows in eq2
    L + mesh.size, ...
    ... c-rows in eq3
    L + 2*mesh.size, ...
    ... y-rows in eq3
    L + 2*mesh.size, ...
    ... G-rows in eq4
    L + 3*mesh.size, ...
    ... x-rows in eq4
    L + 3*mesh.size, ...
    ... y-rows in eq5
    L + 4*mesh.size, ...
    ... x-rows in eq5
    L + 4*mesh.size];
%cols
j = [... p-cols in eq1
    L1,L2,L3,L4, ...
    L1_p, L2_n, L3_p, L4_n, ...
    ... c-cols in eq1
    [L1,L2,L3,L4, ...
    L1_p, L2_n, L3_p, L4_n] + mesh.size, ...
    ... c-cols in eq2
    [L1,L2,L3,L4, ...
    L1_p, L2_n, L3_p, L4_n] + mesh.size, ...
    ... G-cols in eq2
    L + 3*mesh.size, ...
    ... c-cols in eq3
    L + 1*mesh.size, ...
    ... y-cols in eq3
    L + 2*mesh.size, ...
    ... G-cols in eq4
    L + 3*mesh.size, ...
    ... x-cols in eq4
    L + 4*mesh.size, ...
    ... y-cols in eq4
    L + 2*mesh.size, ...
    ... x-cols in eq4
    L + 4*mesh.size];
%vals
v = [... p-coefs in eq1
    -S_r(:);-S_r(:);-SBB_z(:);-SBB_z(:); ...
    S_r(:);S_r(:);SBB_z(:);SBB_z(:); ...
    ... c-coefs in eq1
    -P_r(:);-P_r(:);-PBB_z(:);-PBB_z(:); ...
    P_r(:);P_r(:);PBB_z(:);PBB_z(:); ...
    ... c-coefs in eq2
    q_z_pos(:);q_z_neg(:);-q_z_neg(:);-q_z_pos(:); ...
    qB_r_pos(:);qB_r_neg(:);-qB_r_neg(:);-qB_r_pos(:); ...
    ... G-coefs in eq2
    -mesh.V(:)/dt; ...
    ... c-coefs in eq3
    ones(size(mesh.V(:))); ...
    ... y-coefs in eq3
    ones(size(mesh.V(:)))/dt; ...
    ... G-coefs in eq4
    ones(size(mesh.V(:))); ...
    ... x-coefs in eq4
    -dGdx(x, params); ...
    ... y-coefs in eq5
    ones(size(mesh.V(:))); ...
    ... x-coefs in eq5
    -x(:) ...
    ]';

A = sparse(i,j,v,m,m,numel(v));
%% BC in matrix
% c_in = 0
L = linear_index(J(1), I, mesh);
shift_c = mesh.size;
for l = L(2:end-1)
    A(l+shift_c, l+shift_c) = 1;
    A(l+shift_c, l+3*shift_c) = 0;
    A(l+shift_c, l+1+shift_c) = 0;
    
    A(l+shift_c, l-Nz+shift_c) = 0;
    A(l+shift_c, l+Nz+shift_c) = 0;
end
l = L(1);
    A(l+shift_c, l+shift_c) = 1;
    A(l+shift_c, l+3*shift_c) = 0;
    A(l+shift_c, l+1+shift_c) = 0;
    
    A(l+shift_c, l+Nz+shift_c) = 0;

l = L(end);
    A(l+shift_c, l+shift_c) = 1;
    A(l+shift_c, l+3*shift_c) = 0;
    A(l+shift_c, l+1+shift_c) = 0;
    
    A(l+shift_c, l-Nz+shift_c) = 0;
% p_out = 0
L = linear_index(J(end), I, mesh);
shift_c = mesh.size;
for l = L(2:end-1)
    A(l, l) = 1;
    A(l, l+shift_c) = 0;
    A(l, l-1) = 0;
    A(l, l-1+shift_c) = 0;
    
    A(l, l-Nz) = 0;
    A(l, l-Nz+shift_c) = 0;
    A(l, l+Nz) = 0;
    A(l, l+Nz+shift_c) = 0;
end
l = L(1);
    A(l, l) = 1;
    A(l, l+shift_c) = 0;
    A(l, l-1) = 0;
    A(l, l-1+shift_c) = 0;
    
    A(l, l+Nz) = 0;
    A(l, l+Nz+shift_c) = 0;

l = L(end);
    A(l, l) = 1;
    A(l, l+shift_c) = 0;
    A(l, l-1) = 0;
    A(l, l-1+shift_c) = 0;
    
    A(l, l-Nz) = 0;
    A(l, l-Nz+shift_c) = 0;

d_s = A\rhs;
end

