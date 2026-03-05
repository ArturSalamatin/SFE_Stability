function [x,xi, Psi, X, Phi, Gamma] = calc_full_inlet_solution_asymptotics(...
    x_left, params, sigma, nterms)

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

if(nargin == 3)
    nterms = 10;
end

    Psi = 0;
    X = 0;
    Phi = 0;
    Gamma = 0;
    
if(nterms >= 0)
    Psi = psi0;
    X = b0;
    Phi = phi0;
    Gamma = gamma0;
end
if(nterms >= 1)
    Psi = Psi + t.*psi1;
    X = X + t.*b1;
    Phi = Phi + t.*phi1;
    Gamma = Gamma + t.*gamma1;
end
if(nterms >= 2)
    Psi = Psi + (t.^2).*psi2;
    X = X + (t.^2).*b2;
    Phi = Phi + (t.^2).*phi2;
    Gamma = Gamma + (t.^2).*gamma2;
end
if(nterms >= 3)
    Phi = Phi + (t.^3).*phi3;
    Gamma = Gamma + (t.^3).*gamma3;
end
if(nterms >= 4)
    Phi = Phi + (t.^4).*phi4;
end

Psi = t.^r.*Psi;
X = t.^r.*X;
Phi = t.^r.*Phi;
Gamma = t.^r.*Gamma;

x = x_left;
xi = z_of_x(x,params)/params.z2;
end