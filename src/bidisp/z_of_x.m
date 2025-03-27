function out = z_of_x(x, params)

t = params.t;
a0 = params.a0;
g0 = params.g0;
g1 = params.g1;
r = params.r;

a = params.a;

if(a < a0)
    % dust particles are not extracted yet
    % at the inlet cross-section
    out = (a-x)/g0;
else
    % dust particles are depleted
    % at the inlet cross-section
    out = zeros(size(x));
    % and two sections within the extraction zone exist
    mask = x <= a0;
    out(mask) = (a0 - x(mask))/g0 + z0(params);
        
    mask = (x > a0) & (x <= 1);
    out(mask) = (a-x(mask)+...
        r/(1-r)*log((r+(1-r)*x(mask))/(r+(1-r)*a)))/(1-r);
    
    mask = x > 1;
    if(any(mask))
        error('X must not be greater than 1');
    end
end