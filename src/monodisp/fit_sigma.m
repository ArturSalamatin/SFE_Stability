function [sigma, sol] = fit_sigma(starter, guess, params)
sigma = guess;
[sigma, val] = ...
    fzero(@(s) func_to_min(starter, s, params), guess);
if(abs(val) > 1e-10)
    sigma = NaN;
%     error('Wrong value found!');
end
[val]
if(nargout == 2)
    sol = starter(sigma, params);
    plot_solution(sol.t,sol.y,params.pen)
end
end

function [out, sol] = func_to_min(starter, sigma, params)
sol = starter(sigma, params);
out = sol.y(end,5);
end



