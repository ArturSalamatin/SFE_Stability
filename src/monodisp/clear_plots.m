function clear_plots()
global fig_id

for i = 1:6
    figure(fig_id+i)
    %     legend
    box on
    hold off
end
end