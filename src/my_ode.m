function dy = my_ode(z, y, params, sigma)
% y = (PHI; PSI; GAMMA)

dy = Jac(z, [], params, sigma)*y;

end





