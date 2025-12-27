function [x,xi, Phi, Z, r, psi0, t] = tr_calc_inlet_solution_assymptotics(...
    x_left, params, sigma)
%% assymptotics for [Psi, X] at the jump point, xi = xi0-0
a0 = params.a0;
% x = a0;% linspace(a0, 2*a0, 101);
% xi = z_of_x(x, params)/sol.base_state.z2;
g1 = params.g1;
a = params.a;
alpha = params.r;
C2 = params.C2;

r = (2+sigma)/C2;

psi0 = 0;
b0 = 1;

D = alpha+(1-alpha)*a;
y2 = g1;
m0 = a*D;
m1 = -(D+y2*a);
z1 = a/D;
z2 = -alpha/(2*y2*y2)*(y2/D)^2;
l1 = m0*z1;
l2 = m0*z2+m1*z1;

psi1 = -b0*g1/(D*(r+1));
b1 = (a*(a*psi1-psi0) - b0*(a*a*y2/D+2*a*(1+sigma))-C2*l2*b0*r)/...
    (C2*l1*(r+1)-a*a*(2+sigma));

m2 = y2;
z3 = -alpha/(3*y2*y2)*(y2/D)^3;
l3 = m0*z3+m1*z2+m2*z1;

psi2 = -(a*g1/D*(b1*D-b0*y2)+g1*(a*psi1-psi0)-(D+y2*a)*(r+1)*psi1)/...
    (a*D*(r+2));
b2 = (b0*(1+sigma)+a*(a*psi2-psi1)-b1*(a*a*y2/D+2*a*(1+sigma)) ...
    -C2*(l2*b1*(r+1)+l3*b0*r))/...
    (C2*(r+2)-a*a*(2+sigma));

t = a-x_left;% linspace(0,3e-1,1001);

phi = psi1+psi2*t;
z = b0 + (b1+b2*t).*t;

x = x_left;
xi = z_of_x(x, params)/params.z2;
Phi = phi;
Z = z;
end