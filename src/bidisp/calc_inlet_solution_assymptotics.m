function [x,xi, Psi, X] = calc_inlet_solution_assymptotics(...
    x_left, params, sigma)

[x,xi, Phi, Z, r, psi0, t] = tr_calc_inlet_solution_assymptotics(...
    x_left, params, sigma);

Psi = (t.^r).*(psi0 + t.*Phi);
X = (t.^r).*Z;
end