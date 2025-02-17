function out = Jac_plus(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;

out = zeros(4,4);

c = exp(-y);

out(1,4) = c/2;

out(2,1) = -c/2;
out(2,2) = -1;
out(2,3) = -1;

b = 1/(2-c);
out(3,2) = 2*b;
out(3,3) = (2+(1+sigma)*c)*b;

b = c/2;
out(4,1) = a*a*alpha2*b;
out(4,2) = out(4,1)*R;
out(4,4) = -R*b;

end
