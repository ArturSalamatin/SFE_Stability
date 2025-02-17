function out = Jac_exp(y, ~, sigma, params)

R = params.R;
f = params.alpha*params.a;


out = zeros(4,4);

ksi = exp(y)/(1+exp(y));

out(1,4) = ksi*(1-ksi);

out(2,1) = -ksi*(1-ksi);
out(2,2) = -ksi;
out(2,3) = -ksi;

out(3,2) = 1;
out(3,3) = (1-ksi)*(exp(y) + 2 + sigma);

out(4,1) = ksi*(1-ksi)*f^2;
out(4,2) = out(4,1)*R;
out(4,4) = -R*ksi*(1-ksi);


end