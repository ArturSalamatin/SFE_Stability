function out = z0(params)

a = params.a;
a0 = params.a0;
a1 = params.a1;
r = params.r;
g1 = params.g1;

out = zeros(size(a));
% here z0 does not propagate yet
% out = zeros(size(a));

% here z0 > 0, and popagates
mask = a >= a0;
out(mask) = (a(mask)-a0+...
    r/g1*log((r+g1*a0)./(r+g1*a(mask))))/g1;
end

