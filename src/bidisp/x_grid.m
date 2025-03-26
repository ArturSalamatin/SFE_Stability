function out = x_grid(params, step)
t = params.t;
a0 = params.a0;
a1 = params.a1;
a = params.a; % = sqrt(2t)

if(t<=0)
    error('Time must be positive!');
end
if(a0 > a1)
    error('Dust particle size must be smaller than 1!');
end

if(nargin == 1)
    step = a0/15;
end

if(a <= a0)
    % dust particles are not extracted yet
    % at the inlet cross-section
    n = max(10, ceil(a/step));
    % create a uniform mesh on [0; a] segment
    out = linspace(0, a, n);
else
    n = max(10, ceil(a0/step));
    out = linspace(0, a0, n);
    if(a <= a1)
        % dust particles are depleeted,
        % while large particles --- are not
        % at the inlet cross-section
        step = a0/n;
        n = ceil((a-a0)/step);
        % concatanate meshes
        out = [out(1:end-1), linspace(a0, a, n)];
    else
        error('Large times, t > 0.5, are not considered.')
        step = (1-a0)/n;
        n = ceil((1-a0)/step);
        out = [out(1:end-1), linspace(a0, 1, n)];
        
        n = ceil((a-1)/step);
        out = [out(1:end-1), linspace(1, a, n)];
    end
end




