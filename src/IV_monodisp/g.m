function out = g(x)

if(any(x < 0))
    error('X must be positive!');
end

out = zeros(size(x));
mask = x <= 1;
out(mask) = 1;
end