function out = x_of_z(z, t, params)
% inverse function for z(x)
z0 = params.z0;
a = sqrt(2*t);
a0 = params.a0;
g0 = params.g0;
g1 = params.g1;
r = params.r;

if(a < a0)
    % at initial time frame the dependence
    % is linear. The inversion is simple
    out = a - g0*z;
else
    out = zeros(size(z));
    
    for i = 1:numel(z)
        if(z(i) > z0)
            % a > a0
            % dependence is also linear
            out(i) = a0 - g0*z(i) + ...
                g0/g1*((a - a0 + r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*a))));
        elseif(z == 0)
            out(i) = a;
        else % 0 < z < z0
            % non-linear dependence
            I = find(params.z > z(i), 1, 'last');
            x_guess = (params.x(I) + params.x(I+1))/2;
            out(i) = fzero(@(X) fun(X, a, z(i), r), x_guess);
        end
    end
end
end

function out = fun(X, a, z, r)
out = a - X + r/(1-r)*log((r+(1-r)*X)/(r+(1-r)*a)) - z*(1-r);
end