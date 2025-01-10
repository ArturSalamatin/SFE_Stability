function out = x_grid(params, step)
t = params.t;
a = params.a;

if(t <= 0)
    error('Time must be positive!');
end

if(nargin == 1)
    step = a/105;
end

n = max(10, ceil(a/step));
out = linspace(0, a, n);


end



