function out = X_ass_0(t,z,r)

%% z/z2 -> 0
Z2 = z2(t,r);
xi = z/Z2;
x0 = 1;
x1 = -Z2/(2*t)*(r+(1-r)*sqrt(2*t));
x2 = -(x1.^2)*(r/2)/(r+(1-r)*sqrt(2*t));
x3 = -x1.*x2.*r/(r+(1-r)*sqrt(2*t)) + ...
    x1.^3*(r*(1-r)*sqrt(2*t)/3)*(r+(1-r)*sqrt(2*t)).^(-2);

out = sqrt(2*t)*(x0+xi.*(x1+xi.*(x2+xi.*x3)));


end