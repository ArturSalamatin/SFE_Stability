function [out, dout] = init_value(z, omega)
% perturbation amplitude. Scales the perturbation,
% so it can be set to unity
A = 1;
% cosine perturbation
chi = A*(1-cos(omega*z));

out = [zeros(size(z)), (chi.^2)/2]';
dout = [zeros(size(z)), A*omega*chi.*sin(omega*z)]';

end