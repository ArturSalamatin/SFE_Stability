function functional_plot()
clc

N = 2001;
mesh = uniform_mesh(0, 1, N);
starter = @(sigma, params) starter_omega_Y(sigma, params, mesh);

sigma = linspace(-1.93,0,101);
params.R = 0.2;

f = [350];

for j = 1:numel(f)
    params.f = f(j);    
    out = zeros(size(sigma));
    
    for i = 1:numel(sigma)
        sol = starter(sigma(i), params);
        out(i) = sol.y(end, 5);
        
        %     pen = set_pen('r', '-');
        %     plot_solution(sol.t,sol.y,pen);
    end
    
    figure(9)
    hold on
    box on
    axis([-Inf Inf -Inf Inf])
    xlabel('{\sigma}')
    ylabel('{\Omega(1)}')
    % title()
    
    plot(sigma, out, 'r--', 'LineWidth', 1)
end

end