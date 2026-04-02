function [sigma, sol] = fit_sigma(starter, params, guess)
if(nargin == 2)
guess = make_guess(starter, params);
end
sigma = guess;
[sigma, val] = ...
    fzero(@(s) func_to_min(starter, s, params), guess);
if(abs(val) > 2e-8)
    sigma = NaN;
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

function out_I = make_guess(starter, params)
sigma = linspace(-2,2,101);
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

out_I = (sigma(I) + sigma(I-1))/2;
end
