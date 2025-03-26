function functional_plot()
% the method plots Omega(1, sigma) at various values of R and f
%%
clc
% close all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)
%%
N = 401;
mesh = uniform_mesh(0, 1, N);
starter = @(sigma, params) starter_omega_Y(sigma, params, mesh);

sigma = linspace(-2,0,201);
R = [0.5, 1.3];% [0.4, 0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5];%[0.1,0.2,0.3,0.5,0.7,1.0];

params.f = 5;% 80;% [0.1, 0.3, 0.8, 1];%, 2, 5, 15];

for j = 1:numel(R)
    params.R = R(j);    
    out = zeros(size(sigma));
    
    for i = 1:numel(sigma)
        sol = starter(sigma(i), params);
        out(i) = sol.y(end, 5);
        
        %     pen = set_pen('r', '-');
        %     plot_solution(sol.t,sol.y,pen);
    end
    
    idx = right_monotone(out);
    
    figure(10)
    hold on
    box on
    grid on
    axis([-2 -1 -0.01 0.01])
    xlabel('{\sigma}')
    ylabel('{\Omega(1)}')
    % title()
    
    plot(sigma(idx:end), out(idx:end), 'k--', 'LineWidth', 1)
end

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