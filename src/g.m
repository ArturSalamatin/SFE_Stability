function out = g(x)

if(x < 0)
    error('X must be positive!');
end

if(x <= 1)
    out = 1;
else
    out = 0;
end
end