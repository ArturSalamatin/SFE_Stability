function sol = series_expansion(t, K, sigma, params)
%[Phi, Omega, X, Gamma]

if(sigma == -2)
    error('sigma = -2, expansion is undefined!')
end
if(sigma == -1)
    error('sigma = -1, expansion is undefined!')
end

f2 = (params.f)^2;
R = params.R;

phi0 = 0;
omega0 = 0;
x0 = omega0/(2+sigma);
gamma0 = 1;
c0 = [phi0,omega0,x0,gamma0];

phi1 = gamma0;
omega1 = -phi0 +(2+sigma)*x0;
x1 = -(omega1 - (2+sigma)*x0)/(1+sigma);
gamma1 = f2*phi0 + f2*R*omega0-R*gamma0;

c1 = [phi1,omega1,x1,gamma1];

out = c0 + c1.*t;

T = t;
for i = 2:K
    c2 = coefs(c0,c1,i,params, sigma);
    T = T.*t;
    out = out + c2.*T;
    
    c0 = c1;
    c1 = c2;    
end

sol.y = out;
sol.x = t;

sol = transform(sol);
end

function out = coefs(c0, c1, k, params, sigma)
%[Phi, Omega, X, Gamma]
f2 = (params.f)^2;
R = params.R;

out = zeros(size(c0));

out(1) = c1(4)/k;
out(2) = ((2+sigma)*c1(3)-c1(1))/k;
out(3) = (out(2)+(k-3-sigma)*c1(3))/(k-2-sigma);
out(4) = ((f2*(c1(1)+R*(c1(2)-c0(3))) - R*c1(4)))/k;

end

function sol = transform(sol)
t = sol.x;
% in
%[Phi, Omega, X, Gamma]
Phi = sol.y(:,1);
Omega = sol.y(:,2);
X = sol.y(:,3);
Gamma = sol.y(:,4);
% out
%[Phi, Psi, X, Gamma, Omega, Y, Psi+X]
Y = t.*X;
Psi = Omega - Y;

sol.y = [Phi, Psi, X, Gamma, Omega, Y, Psi+X];
end


