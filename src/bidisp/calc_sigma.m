function [out, sol] = calc_sigma(h, R, starter, pens, params)
global sigma_fig
out = zeros(numel(h), numel(R)) - 1.888;
for j = 1:numel(h)
    j
    params.h = h(j);    
    params.R = R(1);
    params.pen = set_pen(...
        pens.lc{min(1, numel(pens.lc))}, ...
        pens.style{min(j, numel(pens.style))});
    
        if(nargout == 2)
            [sigma, sol] = fit_sigma(starter, params, -params.C1-params.C2);        
        else
            sigma = fit_sigma(starter, params, -params.C1-params.C2);
        end
    out(j,1) = sigma;
    for i = 2:numel(R)
        params.R = R(i);
        I = max(1,i-1);
        
        params.pen = set_pen(...
            pens.lc{min(i, numel(pens.lc))}, ...
            pens.style{1});
        
        if(nargout == 2)
            [sigma, sol] = fit_sigma(starter, params, out(j,I));        
        else
            sigma = fit_sigma(starter, params, out(j,I));
        end
        out(j,i) = sigma;
    end
    
    figure(sigma_fig)
    hold on
    box on
    xlabel('{\itR}')
    ylabel('{\it\sigma}')
    plot(R, out(j,:), 'k-', 'LineWidth', 1)
    
    %         save_figs(sigma, params);
    
end
    figure(sigma_fig)
    plot(0, out(1,1), 'k', ...
        'LineStyle', 'none', ...
        'Marker', params.marker, ...
        'MarkerFacecolor', 'black')
end