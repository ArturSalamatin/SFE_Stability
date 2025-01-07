function dy = my_ode(z, ~, params, sigma)
% y = (PHI; PSI; GAMMA)

dy = Jac(z, [], params, sigma)*y;

end





