function out = G2X(x)
% The function calculates the ratio G(X)/X

out = zeros(size(x));
for i = 1:numel(x)
    if(x(i) < 0)
        error('X must not be negative!');
    elseif(x(i) <= 1)
        out(i) = 1;
    else
        out(i) = 1/x(i);
    end    
end
end