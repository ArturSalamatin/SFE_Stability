function out = z_of_x(x, params)

t = params.t;
a0 = params.a0;
g0 = params.g0;
g1 = 1-params.r;
r = params.r;

a = sqrt(2*t);

if(a < a0)
    out = (a-x)/g0;
else
    out = zeros(size(x));
    mask = x <= a0;
    out(mask) = (a0 - x(mask))/g0 + ...
        (a - a0 + r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*a)))/g1;
    
    mask = (x > a0) & (x <= 1);
    out(mask) = (a-x(mask)+...
        r/(1-r)*log((r+(1-r)*x(mask))/(r+(1-r)*a)))/(1-r);
    
    mask = x > 1;
    if(any(mask))
        error('X must not be greater than 1');
    end
end