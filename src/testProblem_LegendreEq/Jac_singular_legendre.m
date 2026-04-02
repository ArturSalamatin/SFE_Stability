function out = Jac_singular_legendre(x,~,n)

out = zeros(2,2);

lambda = n*(n+1);

out(1,2) = 1;
out(2,1) = -lambda/(1-x^2);
out(2,2) = 2*x/(1-x^2);

end