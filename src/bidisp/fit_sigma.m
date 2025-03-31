function [sigma, sol] = fit_sigma(starter, params, guess)
if(nargin == 2)
    % make guess, if no guess provided for sigma
    [guess, L, R] = make_guess(starter, params);
end
sigma = guess;
[sigma, val, exitflag, output] = ...
    fzero(@(s) func_to_min(starter, s, params), guess);
if(abs(val) > 2e-5)
    % solution may not be found
    sigma = NaN;
    if(nargin == 3)
        % if no guess was constructed internally
        [guess, L, R] = make_guess(starter, params);
    end
    % use segment division by half method
    [sigma, val] = ...
        fzero(@(s) func_to_min(starter, s, params), [L,R]);
    %     error('Wrong value found!');
end

if(nargout == 2)
    sol = starter(sigma, params);
    plot_solution(sol.t,sol.y,params.pen)
end
end

function [out, sol] = func_to_min(starter, sigma, params)
sol = starter(sigma, params);
out = sol.y(end,5);
end

%%
function out = right_monotone(y)
% identify monotone interval in y-values
i = numel(y);
while (i > 1) && (y(i) > y(i-1))
    i = i-1;
end
out = i;
end

function out = sign_change(y)
% identify the change of sign in y-values
i = 1;
I = numel(y);
while (i <= I) && (y(i) < 0)
    i = i+1;
end
out = i;
end

function [out_I, L, R] = make_guess(starter, params)
a = params.a;
sigma_min = (-1-a*a*params.dz2dt/params.z2)*1.1;

sigma = linspace(sigma_min,1,81);
out = zeros(size(sigma));
for i = 1:numel(sigma)
    sol = starter(sigma(i), params);
    out(i) = sol.y(end, 5);
end
idx = right_monotone(out);
% identify interval of monotonicity
sigma = sigma(idx:end);
out = out(idx:end);

I = sign_change(out);
R = sigma(I);
L = sigma(I-1);
out_I = (L + R)/2;
end
