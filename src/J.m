function out = J(x, params)

out = 0;

R = params.R;
alpha2 = (params.alpha)^2;
a2 = (params.a)^2;

sigma = x(1);
phi_plus = x(2);
chi_plus = x(3);
g_plus = x(4);

bb = -1;
IC = [0;0;0;bb];
A = -30;

c1 = exp(A);
phi0 = 0.5;
psi0 = 0;
chi0 = 0;
gamma0 = -R/2;

c2 = c1*c1;
phi1 = -R/8;
psi1 = 1/8;
chi1 = -1/(8*sigma);
gamma1 = -a2*alpha2/8*(1+R*R);

IC = IC + ...
    bb*[phi0; psi0; chi0; gamma0]*c1 + ...
    bb*[phi1; psi1; chi1; gamma1]*c2;

Jac = @(z, u) Jac_minus(z, u, params, sigma);
[z_minus, u_minus] = solver(sigma, params, Jac, IC, A);

% IC = [phi_plus; 0; chi_plus; g_plus];
% Jac = @(z, u) Jac_plus(z, u, params, sigma);
% [z_plus, u_plus] = solver(sigma, params, Jac, IC, -A);

% out = ...log10...
%     ( ...
%     ... abs(z2 - z_end) + ...
%     sum(abs(u_minus(end, :)-u_plus(end,:)).^2));

% z = [exp(z_minus)/2; 1-exp(-z_plus(end:-1:1))/2];
% u = [u_minus;        u_plus(end:-1:1, :)];

end



function out = Jac_plus(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;

out = zeros(4,4);

c = exp(-y);

out(1,4) = c/2;

out(2,1) = -c/2;
out(2,2) = -1;
out(2,3) = -1;

b = 1/(2-c);
out(3,2) = 2*b;
out(3,3) = (2+(1-sigma)*c)*b;

b = c/2;
out(4,1) = a*a*alpha2*b;
out(4,2) = out(4,1)*R;
out(4,4) = -R*b;

end




