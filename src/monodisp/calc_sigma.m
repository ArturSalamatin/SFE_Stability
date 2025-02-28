function out = calc_sigma(f, R, starter, pens)
out = zeros(numel(f), numel(R)) - 1.888;
for j = 1:numel(f)
    params.f = f(j);
    for i = 1:numel(R)
        params.R = R(i);
        I = max(1,i-1);
        
        pen = set_pen(pens.lc{min(i, numel(pens.lc))}, pens.style);
        params.pen = pen;
        
        [sigma, ~] = fit_sigma(starter, out(j,I), params);
        out(j,i) = sigma;
        
        
        
        
        %         x = sol.t';
        %         x = x(x<0.2);
        %         sol = series_expansion(x, 805, sigma, params);
        %
        %         pen_series.lc = 'm';
        %         pen_series.style = '--';
        %         plot_solution(sol.x,sol.y,pen_series);
    end
    
    figure(7)
    hold on
    box on
    xlabel('{\itR}')
    ylabel('{\it\sigma}')
    plot(R, out(j,:), 'k-', 'LineWidth', 1)
    
    %         save_figs(sigma, params);
    
end


end