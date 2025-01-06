function dy = my_ode(xi, y, t, r, R, alpha, sigma, xi_max)
% y = (PHI; PSI; GAMMA)

dy = Jac(xi, y, t, r, R, alpha, sigma, xi_max)*y;

end





