function [x,xi, Psi, X] = calc_solution_assymptotics(params, sol)
%% assymptotics for [Psi, X] at the jump point, xi = xi0-0
a0 = params.a0;
% x = a0;% linspace(a0, 2*a0, 101);
% xi = z_of_x(x, params)/sol.base_state.z2;
xi0 = z_of_x(a0, params)/sol.base_state.z2;
g1 = params.g1;
g0 = params.g0 - g1;
a = params.a;
alpha = params.r;
dz0dt = sol.base_state.dz0dt;
C1 = sol.base_state.C1;
C2 = sol.base_state.C2;
sigma = sol.sigma;


psi0 = g0*a*dz0dt - (g0+g1)/C1*(C2 + (1+sigma)*(1-xi0));
b0 = 1;

D = alpha+(1-alpha)*a;
y3 = alpha+(1-alpha)*a0;
y0 = a0*y3;
y1 = 2*y3-alpha;
y2 = 1-alpha;

z0 = (a-a0)/y2+alpha/(y2^2)*log(y3/D);
l0 = y0*z0;

psi1 = (a*y2*y3/D+a0*y2*psi0)/y0;
b1 = -(a*a/D*y3 + (1+sigma)*a0*a0+a*a0*psi0)/(l0*C2);

z1 = -a0/y3;
l1 = y0*z1+y1*z0;

psi2 = (a*y2/D*(y3*b1+y2)+y2*(a0*psi1+psi0) - psi1*y1)/(2*y0);
b2 = -((a*a/D*(y2+y3*b1)+(1+sigma)*(a0*a0*b1+2*a0)+a*(psi0+a0*psi1))/C2 + b1*l1)/(2*l0);

z2 = -alpha/(2*y3*y3);
l2 = (y2*z0+y1*z1+y0*z2);

psi3 = (a*y2/D*(y3*b2+y2*b1)+y2*(a0*psi2+psi1) - psi1*y2 - 2*psi2*y1)/(3*y0);
b3 = -((a*(psi2*a0+psi1) + a*a/D*(y2*b1+y3*b2) + (1+sigma)*(1+2*a0*b1+a0*a0*b2))/C2+b1*l2+2*b2*l1)/(3*l0);

z3 = alpha*y2/(3*y3*y3*y3);
l3 = z3*y0 + z2*y1 + z1*y2;
psi4 = (a*y2/D*(b3*y3+b2*y2) + y2*(a0*psi3+psi2) - 3*psi3*y1-2*psi2*y2)/(4*y0);
b4 =  -(...
    (...
    a*a/D*(b3*y3+b2*y2)...
    +(1+sigma)*(b3*a0*a0+b2*2*a0+b1)...
    + a*(psi3*a0+psi2)...
    )/C2 ...
    + l3*b1+2*b2*l2+3*b3*l1)/(4*l0);

u = linspace(0,a0,101);

psi = psi0+(psi1+(psi2+(psi3+psi4*u).*u).*u).*u;
b = b0 + (b1+(b2+(b3+b4*u).*u).*u).*u;

x = a0+u;
xi = z_of_x(x, params)/sol.base_state.z2;
Psi = psi;
X = b;
end