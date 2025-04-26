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

% figure(3000)
% hold on
% plot(sigma, 0, 'sk', 'MarkerFaceColor', 'black')

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
flag = (y(i) > y(i-1));
while (i > 1) && ((y(i) > y(i-1)) == flag)
    i = i-1;
end
out = i;
end

function out = sign_change(y)
% identify the change of sign in y-values
i = 1;
I = numel(y);
while (i <= I-1) && (y(i)*y(i+1) > 0)
    i = i+1;
end
out = i;
end

function [out_I, L, R] = make_guess(starter, params)
sigma_min = -params.C1-params.C2;

sigma = linspace(-4,-0.6,31);
out = zeros(size(sigma));
for i = 1:numel(sigma)
    sol = starter(sigma(i), params);
    out(i) = sol.y(end, 5);
end

% figure(3000)
% hold on
% % axis([-Inf Inf -1 1])
% plot(sigma, out, 'k-', 'LineWidth', 1)
% hold on


idx = right_monotone(out);
% figure(3000)
% plot(sigma(idx:end), out(idx:end), 'r-', 'LineWidth', 1)
% identify interval of monotonicity
out_R = out(idx);

if(out_R*out(end) < 0)
    % a nice segment of sign change is found
    sigma = sigma(idx:end);
    out = out(idx:end);
    I = sign_change(out);
    R = sigma(I);
    L = sigma(I+1);
    out_I = (L + R)/2;
    return
end

% otherwise, better approx near singular point required
I = [idx-1, idx];
[s, o] = get_sign_change(...
    sigma(I), out(I), starter, params);
    sigma = sigma(idx:end);
    out = out(idx:end);
sigma = [s,sigma];
out = [o,out];

% figure(3000)
% hold on
% axis([-Inf Inf -5 5])
% plot(sigma, out, 'k-', 'LineWidth', 1)
% plot(sigma_min, 0, 'ok', 'MarkerFaceColor', 'black')

I = sign_change(out);

R = sigma(I);
L = sigma(I+1);
out_I = (L + R)/2;
end

function [s, o] = get_sign_change(sigma, out, starter, params)

l = sigma(1);
r = sigma(2);
val_l = out(1);
val_r = out(2);

while(val_l*val_r > 0)
    c = (l+r)/2;
    sol = starter(c, params);
    val = sol.y(end, 5);
    if(val < 0)
        s = c;
        o = val;
        return
    end
        
    if(val > val_l)
        val_l = val;
        l = c;
    else
        val_r = val;
        r = c;
    end
end
s = l;
o = val_l;
end