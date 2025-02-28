function plot_solution(t,y,pen)
fig_id = 700;

factor = 1;

names = {'{\Phi}','{\Psi}','X','{\Gamma}','{\Omega}','Y', '{\Psi}+X'};

for i = 1:7
    figure(fig_id+i)
    box on
    %     axis([0 1 -Inf Inf])
    plot(t,y(:,i)/factor...
        , 'LineWidth', 1 ...
        , 'Color', pen.lc ...
        , 'LineStyle', pen.style ...
        )
    hold on
%     plot(t([1,end]), y([1,end],i)/factor, 'o', 'MarkerFaceColor', 'black')
    xlabel('{\xi}')
    ylabel(names{i})
end
end