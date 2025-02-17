function out = Jac_changed_vars(y, ~, params, sigma)

if(nargin == 3)
    sigma = params.sigma;
end


a = params.a;
R = params.R;
alpha2 = (params.alpha)^2;

f2 = a*a*alpha2;

out = zeros(4,4);

out(1,4) = 1;

out(2,1) = -1;
out(2,3) = 2+sigma;

out(3,2) = 1/(y*(1-y));
out(3,3) = (2+sigma)/y;

out(4,1) = f2;
out(4,2) = out(4,1)*R;
out(4,3) = -out(4,2)*y;
out(4,4) = -R;

end