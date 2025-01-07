function out = x_of_z(z, params)

z0 = params.z0;
a = params.a;
a0 = params.a0;
g0 = params.g0;
g1 = params.g1;
r = params.r;

if(a < a0)
    out = a - g0*z;
else
    out = zeros(size(z));
    
    for i = 1:numel(z)
        if(z(i) > z0)
            % a > a0
            out(i) = a0 - g0*z(i) + ...
                g0/g1*((a - a0 + r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*a))));
        elseif(z == 0)
            out(i) = a;
        else % 0 < z < z0
            I = find(params.z > z(i), 1, 'last');
            x_guess = (params.x(I) + params.x(I+1))/2;
            out(i) = fzero(@(X) fun(X, a, z(i), r), x_guess);
        end
    end
end
end

function out = fun(X, a, z, r)
out = a - X + r/(1-r)*log((r+(1-r)*X)/(r+(1-r)*a)) - z*(1-r);
end

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