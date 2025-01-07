function out = Jac(z, y, sigma, params)

a = params.a;
R = params.R;
alpha = params.alpha;


X0 = x0(z, params);


dxdt = 1/G(a)*G2X(X0, params);
dcdz = -g(X0, params)*dxdt;

out = zeros(3,3);
out(1,3) = 1;
out(2,1) = dcdz;
out(2,2) = -sigma*g(X0, params)./(sigma*X0 + dxdt);
out(3,1) = alpha^2;
out(3,2) = out(3,1)*R;
out(3,3) = - R*dcdz;

end

function out = x0(z, params)
out = 1/(1-r)*(G(sqrt(2*t), r)*exp(-xi*(1-r))-r);
end
