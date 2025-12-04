function [x,xi, Psi, X] = calc_inlet_solution_assymptotics(fig_id,pen, params, sol)
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

r = (2+sigma)/C2;

psi0 = 0;
b0 = 0.76;

D = alpha+(1-alpha)*a;
y2 = g1;
m0 = a*D;
m1 = -(D+y2*a);
z1 = a/D;
z2 = -1/(2*y2)*(y2/D)^2;
l1 = m0*z1;
l2 = m0*z2+m1*z1;

psi1 = -b0*g1/(D*(r+1));
b1 = (a*(a*psi1-psi0) - b0*(a*a*y2/D+2*a*(1+sigma))-C2*l2*b0*r)/...
    (C2*l1*(r+1)-a*a*(2+sigma));

m2 = y2;
z3 = -1/(3*y2)*(y2/D)^3;
l3 = m0*z3+m1*z2+m2*z1;

psi2 = -(a*g1/D*(b1*D-b0*y2)+g1*(a*psi1-psi0)-(D+y2*a)*(r+1)*psi1)/...
    (a*D*(r+2));
b2 = (b0*(1+sigma)+a*(a*psi2-psi1)-b1*(a*a*y2/D+2*a*(1+sigma)) ...
    -C2*(l2*b1*(r+1)+l3*b0*r))/...
    (C2*(r+2)-a*a*(2+sigma));

t = linspace(0,1e-1,101);

psi = (psi0+(psi1+psi2*t).*t).*(t.^r);
b = (b0 + (b1+b2*t).*t).*(t.^r);

x = a-t;
xi = z_of_x(x, params)/sol.base_state.z2;
Psi = psi;
X = b;



figure(fig_id+2)
hold on
plot(xi, Psi ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-.')

figure(fig_id+3)
hold on
plot(xi, X ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-.')
    
    my_figure(fig_id+10)
hold on
    plot(xi, Psi./X ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-.')

end