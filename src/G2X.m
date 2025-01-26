function out = G2X(x)
% The function calculates the ratio G(X)/X


if(any(x < 0))
    error('X must not be negative!');
end

out = ones(size(x));
mask = x > 1;
out(mask) = 1./x(mask);
end