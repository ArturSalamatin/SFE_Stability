function out = X_ass_1(t,z,r)

%% z/z2 -> 1
% d = 1/(r+(1-r)*sqrt(2*t));
Z2 = z2(t,r);
xi = sqrt(1-z/Z2);
x1 = sqrt(r*Z2/t);
x2 = x1*x1*sqrt(2*t)*(1-r)/(6*r);
x3 = (1-r)/r*sqrt(t/2)*(x1.*x2) - (1-r)^2/(4*r*r)*t*x1.^3 - (x2.^2)./(2*x1);
x4 = (1-r)/r*sqrt(2*t)/6*(3*x3.*x1+2*x2.*x3+x2.*x2) + ...
    ((1-r)/r)^3*(2*t).^(5/2)*(x1.^4)/10 - x2.*x3 - ...
    ((1-r)/r)^2*t*x1.^2*x2;

out = sqrt(2*t)*xi.*(x1+xi.*(x2+0*xi.*(x3+xi.*x4)));


end