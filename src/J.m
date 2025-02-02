function out = J(x, params)

sigma = x(1);
phi_plus = x(2);
chi_plus = x(3);
g_plus = x(4);

IC = [0;0;0;1];
A = -20;
Jac = @(z, u) Jac_minus(z, u, params, sigma);
[z_minus, u_minus] = solver(sigma, params, Jac, IC, A);

IC = [phi_plus; 0; chi_plus; g_plus];
Jac = @(z, u) Jac_plus(z, u, params, sigma);
[z_plus, u_plus] = solver(sigma, params, Jac, IC, -A);

out = log10( ...
    ... abs(z2 - z_end) + ...
    sum(abs(u_minus(end, :)-u_plus(end,:))));

z = [exp(z_minus)/2; 1-exp(-z_plus(end:-1:1))/2];
u = [u_minus;        u_plus(end:-1:1)];

end


function out = Jac_minus(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

out = zeros(4,4);

c = exp(y);

out(1,4) = c/2;

b = -c/(4-2*c);
out(2,1) = (2-c)*b;
out(2,2) = b;
out(2,3) = b;

b = 1/(2-c);
out(3,2) = 2*b;
out(3,3) = (6-(1+sigma)*c)*b;

b = c/2;
out(4,1) = a*a*alpha2*b;
out(4,2) = out(4,1)*R;
out(4,4) = -R*b;

end


function out = Jac_plus(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

out = zeros(4,4);

c = exp(-y);

out(1,4) = c/2;

out(2,1) = -c/2;
out(2,2) = -1;
out(2,3) = -1;

b = 1/(2-c);
out(3,2) = 2*b;
out(3,3) = (2-(1-sigma)*c)*b;

b = c/2;
out(4,1) = a*a*alpha2*b;
out(4,2) = out(4,1)*R;
out(4,4) = -R*b;

end




