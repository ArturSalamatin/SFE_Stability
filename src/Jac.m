function out = Jac(z, y, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end


a = params.a;
R = params.R;
alpha = params.alpha;


X0 = x_of_z(z, params);


dxdt = 1/G_of_x(a, params)*G2X(X0, params);
dcdz = -g(X0, params)*dxdt;

out = zeros(3,3);
out(1,3) = 1;
out(2,1) = dcdz;
out(2,2) = -sigma*g(X0, params)./(sigma*X0 + dxdt);
out(3,1) = alpha^2;
out(3,2) = out(3,1)*R;
out(3,3) = - R*dcdz;

end
