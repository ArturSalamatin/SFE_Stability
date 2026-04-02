function out = G_div_X(x, params)
% The function calculates the ratio G(X)/X
r = params.r;
a0 = params.a0;
a1 = params.a1;
g0 = params.g0;
g1 = params.g1;

out = zeros(size(x));
for i = 1:numel(x)
    if(x(i) < 0)
        error('X must not be negative!');
    elseif(x(i) <= a0)
        out(i) = g0;
    elseif(x(i) <= a1)
        out(i) = r/x(i)+g1;
    else
        out(i) = 1/x(i);
    end    
end
end