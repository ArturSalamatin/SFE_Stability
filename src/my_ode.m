function dy = my_ode(z, y, params, sigma, Jac)
% y = (PHI; PSI; CHI; GAMMA)

dy = Jac(z, [])*y;

end





