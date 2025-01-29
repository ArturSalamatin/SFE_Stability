function out = J(sigma, params)
[z, y] = solver(sigma, params);

% z2 = params.z2;

% z_end = z2 - h;
% bc = BC(sigma, psi0, gamma0, z_end, params);

out = log10( ...
    ... abs(z2 - z_end) + ...
    sum(abs(y(end, 2))));
end


function out = BC(sigma, psi0, gamma0, z, params)
% sigma and psi0 are to be fitted
a = params.a;
R = params.R;
alpha = params.alpha;

psi0 = psi0;
phi0 = psi0;
gamma0 = gamma0;

psi1 = -gamma0 - a*sigma*psi0;
phi1 = -gamma0;
gamma1 = (alpha^2*(1+R*sigma)*psi0 - R/a*gamma0);

out = [phi0+phi1*z, (psi0+psi1*z)*z^2, gamma0+gamma1*z];
end




