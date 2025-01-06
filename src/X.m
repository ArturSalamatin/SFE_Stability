function out = X(t, z, r)

x_max = min(X_ass_0(t,z,r), 1);
x_min = max(0, X_ass_1(t,z,r));
out = zeros(size(x_max));
for ii = 1:numel(z)
    if(z(ii) == 0)
        out(ii) = sqrt(2*t);
    elseif(z(ii) == z2(t, r))
        out(ii) = 0;
    else
        out(ii) = fzero(@(X) fun(X, t, z(ii), r), [x_min(ii), x_max(ii)]);   
    end
end

end

function out = fun(X, t, z, r)
out = sqrt(2*t) - X + r/(1-r)*log((r+(1-r)*X)/(r+(1-r)*sqrt(2*t))) - z*(1-r);
end
