function out = x_grid(params, step)
t = params.t;
a0 = params.a0;

if(t<=0)
    error('Time must be positive!');
end
if(a0 > 1)
    error('Dust particle size must be smaller than 1!');
end

if(nargin == 1)
    step = a0/15;
end

a = sqrt(2*t);
if(a < a0)
    n = max(10, ceil(a/step));
    out = linspace(0, a, n);
else
    n = max(10, ceil(a0/step));
    out = linspace(0, a0, n);
    if(a < 1)
        step = a0/n;
        n = ceil((a-a0)/step);
        out = [out(1:end-1), linspace(a0, a, n)];
    else
        step = (1-a0)/n;
        n = ceil((1-a0)/step);
        out = [out(1:end-1), linspace(a0, 1, n)];
        
        n = ceil((a-1)/step);
        out = [out(1:end-1), linspace(1, a, n)];
    end
end




