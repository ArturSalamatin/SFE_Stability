function out = ode_legendre(x,y, jac)

out = jac(x,[])*y;

end