function out = mass_legendre(t,~)

out = eye(2,2);
out(2,2) = 1-t^2;

end

