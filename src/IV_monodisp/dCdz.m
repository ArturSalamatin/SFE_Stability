function dcdz = dCdz(X, t)
dxdt = dXdt(X, t);
dcdz = g(X).*dxdt;
end

function out = dXdt(X, t)
out = G2X(X)/G_of_x(a(t));
end

function out = G2X(x)
% The function calculates the ratio G(X)/X

if(any(x < 0))
        error('X must not be negative!');
end
out = ones(size(x));
mask = x > 1;
out(mask) = 1./x(mask);
mask = x == 0;
out(mask) = 0;
end

function out = G_of_x(x)
out = min(1,x);
end