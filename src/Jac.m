function out = Jac(xi, y, t, r, R, alpha, sigma, xi_max)

% global X0 C0

X0 = x0(xi, t, r);
C0 = c0(X0, t, r);

out = zeros(3,3);
out(1,3) = X0;
out(2,1) = -(1-r)*(1-C0);
out(2,2) = -(1-r)*sigma*X0.^2./(sigma*X0.^2 + 1-C0);
out(3,1) = alpha^2*X0;
out(3,2) = out(3,1)*R;
out(3,3) = - R*(1-r)*(1-C0);

out = out*xi_max;

end

function out = x0(xi, t, r)
out = 1/(1-r)*(G(sqrt(2*t), r)*exp(-xi*(1-r))-r);
end

function out = c0(X0, t, r)
out = 1 - G(X0, r)/G(sqrt(2*t), r);
end