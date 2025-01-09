function out = Jac(z, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

R = params.R;
alpha = params.alpha;

X0 = x_of_z(z, params);

dxdt = dXdt(X0, params);
dcdz = -g(X0, params)*dxdt;

out = zeros(3,3);
out(1,3) = 1;
out(2,1) = -dcdz/sigma;
out(2,2) = (dxdt./X0/sigma-1)*g(X0, params)./X0;
out(3,1) = alpha^2;
out(3,2) = out(3,1)*R*sigma;
out(3,3) = - R*dcdz;

end
