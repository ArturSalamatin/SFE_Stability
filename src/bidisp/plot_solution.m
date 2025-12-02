function plot_solution(t,y,pen)
fig_id = 700;

factor = 1;

names = {'{\Phi}','{\Psi}','X','{\Gamma}','{\Omega}','Y', '{\Psi}+X', 'Q', 'P', 'P^{\prime}'};

for i = [2,3,5,6]% 1:6%numel(names)
    if(i > size(y,2))
        break
    end
    my_figure(fig_id+i)
    box on
    %     axis([0 1 -Inf Inf])
    plot(t,y(:,i)/factor...
        , 'LineWidth', 1 ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-'... pen.style{1}(1) ...
        )
    hold on
%     plot(t([1,end]), y([1,end],i)/factor, 'o', 'MarkerFaceColor', 'black')
    xlabel('{\xi}')
    ylabel(names{i})
end

% figure(fig_id+8)
% axis([-Inf Inf -Inf 5])
end