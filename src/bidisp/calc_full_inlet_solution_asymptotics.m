function [x,xi, Psi, X, Phi, Gamma] = calc_full_inlet_solution_asymptotics(...
    x_left, params, sigma)

a = params.a;
alpha= params.r;
C1 = params.C1;
C2 = params.C2;
zeta2 = params.z2;
g1 = params.g1;
h = params.h;
R = params.R;

D = alpha+(1-alpha)*a;
r = (2+sigma)/C2;

t = a - x_left;

psi0 = 0;
gamma0 = 0;
b0 = 1;
phi0 = 0;

omega0 = a*D;
y2 = 1-alpha;
omega1 = D+a*y2;
m0 = a*D;
m1 = -omega1;
m2 = y2;
z1 = a/D;
z2 = -alpha/(2*D*D);
z3 = -alpha*y2/(3*D*D*D);
l1 = a*a;
l2 = m0*z2+m1*z1;
l3 = m0*z3+m1*z2+m2*z1;

lambda0 = a*a*(2+sigma);
lambda1 = a*a/D*y2+2*a*(1+sigma);
lambda2 = 1+sigma;

phi1 = a*gamma0/(D*zeta2*(r+1));
gamma1 = 0;
psi1 = -g1/(D*(r+1))*(phi0+b0);
b1 = (a*a*psi1 - (C2*l2*r+lambda1)*b0)/(C2*l1*(r+1)-lambda0);

phi2 = (a*gamma1-gamma0+zeta2*y2*(r+1)*phi1)/(zeta2*D*(r+2));
gamma2 = h*h*a/(zeta2*D*(r+2))*(phi1+R*psi1);
psi2 = -(...
    g1/D*(a*D*(phi1+b1)-a*b0*y2-(a+y2)*phi0)+(a*g1-omega1*(r+1))*psi1)/...
    (omega0*(r+2));
b2 = ...
    (a*a*psi2-a*psi1-(lambda1+C2*l2*(r+1))*b1+(lambda2-C2*l3*r)*b0)...
    /(C2*l1*(r+2)-lambda0);

phi3 = a*gamma2/(zeta2*D*(r+3));
gamma3 = ...
    (zeta2*y2*(r+2)*gamma2+...
    h*h*a*(phi2+R*psi2)...
    -h*h*(phi1+R*psi1)-R*g1*a/C1*gamma2)/...
    (zeta2*D*(r+3));

phi4 = (zeta2*y2*(r+3)*phi3+a*gamma3-gamma2)/(zeta2*D*(r+4));


Psi = t.^r.*(psi0 + t.*(psi1+t.*psi2));
X = t.^r.*(b0 + t.*(b1+t.*b2));
Phi = t.^r.*(phi0 + t.*(phi1+t.*(phi2+t.*(phi3+t.*phi4))));
Gamma = t.^r.*(gamma0 + t.*(gamma1+t.*(gamma2+t.*gamma3)));

x = x_left;
xi = z_of_x(x,params)/params.z2;

end