function [x,xi, Psi, X] = calc_inlet_solution_asymptotics(...
    x_left, params, sigma)

[x,xi, Phi, Z, r, psi0, t] = calc_tr_inlet_solution_asymptotics(...
    x_left, params, sigma);

Psi = (t.^r).*(psi0 + t.*Phi);
X = (t.^r).*Z;
end