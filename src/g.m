function out = g(x)

if(x < 0)
    error('X must be positive!');
end
out = x <= 1;
end