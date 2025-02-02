function out = Jac(z, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

R = params.R;
a = params.a;
alpha2 = (params.alpha)^2;

out = zeros(4,4);
out(1,4) = 1;

out(2,1) = -1;
out(2,3) = 1+sigma;

out(3,2) = 1;
out(3,3) = (2+sigma)*(1-z) + z;

out(4,1) = a*a*alpha2;
out(4,2) = out(4,1)*R;
out(4,4) = -R;

end
