function out = G2X(x, a0, r)
% The function calculates the ratio G(X)/X

g0 = 1-r +r/a0;

out = size(x);
for i = 1:numel(x)
    if(x(i) < 0)
        error('X must not be negative!');
    elseif(x(i) < a0)
        out(i) = g0;
    elseif(x(i) < 1)
        out(i) = r/x(i)+1-r;
    else
        out(i) = 1/x(i);
    end    
end
end