function out = Jac(z, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

R = params.R;
alpha2 = (params.alpha)^2;

X0 = x_of_z(z, params);

dxdt = dXdt(X0, params);
dcdz = dCdz(X0, params);

out = zeros(3,3);
out(1,3) = 1;
out(2,1) = -dcdz.*(sigma*X0+dxdt);
out(2,2) = -sigma;
out(3,1) = alpha2;
out(3,2) = out(3,1)*R;
out(3,3) = - R*out(2,1);

end
