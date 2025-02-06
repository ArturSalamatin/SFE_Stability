function out = Jac(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end

a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;

out = zeros(4,4);

out(1,4) = 1;

out(2,1) = -1;
out(2,2) = -1/(1-y);
out(2,3) = -1/(1-y);

out(3,2) = 1/(y*(1-y));
out(3,3) = (2-y+sigma*(1-y))*out(3,2);

out(4,1) = a*a*alpha2;
out(4,2) = out(4,1)*R;
out(4,4) = -R;


end