function out = x_of_z(z, t, params)
% inverse function for z(x)
a0 = params.a0;
g0 = params.g0;
g1 = params.g1;
r = params.r;

tn = numel(t);
zn = numel(z);

[t, z] = meshgrid(t,z);

out = zeros(size(z));

parfor ti = 1:tn
    Z0 = z0(t(1,ti),params);
    Z2 = z2(t(1,ti),params);
    for zj = 1:zn
        A = a(t(zj,ti));
        if(z(zj,ti) > Z2)
            out(zj,ti) = 0;
            continue;
        end
        % z < z2
        if(A < a0)
            % at initial time frame the dependence
            % is linear. The inversion is simple
            out(zj,ti) = A - g0*z(zj,ti);
        else
            % a > a0
            if(z(zj,ti) > Z0)
                % dependence is also linear
                out(zj,ti) = a0 - g0*z(zj,ti) + ...
                    g0/g1*((A - a0 + r/(1-r)*log((r+(1-r)*a0)/(r+(1-r)*A))));
            elseif(z(zj,ti) == 0)
                out(zj,ti) = A;
            else % 0 < z < z0
                % non-linear dependence
                out(zj,ti) = fzero(@(X) fun(X, A, z(zj,ti), r), [a0, A]);
            end
        end
    end
end
end

function out = fun(X, a, z, r)
out = a - X + r/(1-r)*log((r+(1-r)*X)/(r+(1-r)*a)) - z*(1-r);
end