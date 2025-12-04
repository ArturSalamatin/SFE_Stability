function plot_solution_assymptotics(fig_id,pen, params, sol)

[~,xi, Psi, X] = calc_solution_assymptotics(params, sol);

figure(fig_id+2)
hold on
plot(xi, Psi ...
        , 'Color', pen.lc ...
        , 'LineStyle', '--')

figure(fig_id+3)
hold on
plot(xi, X ...
        , 'Color', pen.lc ...
        , 'LineStyle', '--')

end