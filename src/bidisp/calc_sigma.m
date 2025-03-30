function out = calc_sigma(h, R, starter, pens, params)
out = zeros(numel(h), numel(R)) - 1.888;
for j = 1:numel(h)
    params.h = h(j);
    
    params.R = R(1);    
    pen = set_pen(...
        pens.lc{min(1, numel(pens.lc))}, ...
        pens.style{min(j, numel(pens.style))});
    params.pen = pen;
    
    [sigma] = fit_sigma(starter, params);
    out(j,1) = sigma;
    for i = 2:numel(R)
        params.R = R(i);
        I = max(1,i-1);
        
        pen = set_pen(pens.lc{min(i, numel(pens.lc))}, pens.style);
        params.pen = pen;
        
        [sigma] = fit_sigma(starter, params, out(j,I));
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