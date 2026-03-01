function state_s = single_iteration(...
    state_s, state_prev, mesh, params, dt)
%SINGLE_ITERATION Summary of this function goes here
%   Detailed explanation goes here

% ratio of two spatial scales
B = params.B;
% nmbr of unknowns/equations
eq_nmbr = 5;
Nz = mesh.Nz;
Nr = mesh.Nr;

dr = mesh.dr;
dz = mesh.dz;

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
L = linear_index(J_l, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*q_z_pos(:);
% j, i
L = L+1;
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*q_z_neg(:);
%j+1, i
L = linear_index(J_r, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*q_z_neg(:);
%j,i
L = L-1;
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*q_z_pos(:);

%j, i-1
L = linear_index(J, I_l, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*qB_r_pos(:);
%j, i
L = L+Nz;
rhs(L + shift) = ...
    rhs(L + shift) - c(L).*qB_r_neg(:);
%j, i+1
L = linear_index(J, I_r, mesh);
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*qB_r_neg(:);
%j, i
L = L-Nz;
rhs(L + shift) = ...
    rhs(L + shift) + c(L).*qB_r_pos(:);
%% set rhs eqn: div(v) = 0
shift = 0;
% j-1, i
L = linear_index(J_l, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*S_r(:);
% j, i
L = L+1;
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*S_r(:);
%j+1, i
L = linear_index(J_r, I, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*S_r(:);
%j,i
L = L-1;
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*S_r(:);

%j, i-1
L = linear_index(J, I_l, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*SBB_z(:);
%j, i
L = L+Nz;
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*SBB_z(:);
%j, i+1
L = linear_index(J, I_r, mesh);
rhs(L + shift) = ...
    rhs(L + shift) - p(L).*SBB_z(:);
%j, i
L = L-Nz;
rhs(L + shift) = ...
    rhs(L + shift) + p(L).*SBB_z(:);
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
%% problem matrix
m = eq_nmbr*mesh.size;
nnz = 7*mesh.size + 3*5*mesh.size;
% 
% l1 = 
% 
% i = []

A = sparse([],[],[],m,m,nnz);
end

