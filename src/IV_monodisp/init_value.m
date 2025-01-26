function [out, dout] = init_value(z, omega)
% perturbation amplitude. Scales the perturbation,
% so it can be set to unity
A = 1;
% cosine perturbation
chi = A*(1-cos(omega*z));

out = 0*[zeros(size(z)), (chi.^2)/2]';
dout = 0*[zeros(size(z)), A*omega*chi.*sin(omega*z)]';

end