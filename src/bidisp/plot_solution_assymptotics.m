function plot_solution_assymptotics(fig_id,pen, params, sol)
a0 = params.a0;
x_right = a0 + linspace(a0,0,1001); 

[~,xi, Psi, X] = calc_solution_assymptotics(...
x_right, params, sol.sigma);

figure(fig_id+1)
hold on
plot(xi, Psi ...
        , 'Color', pen.lc ...
        , 'LineStyle', '--')

figure(fig_id+2)
hold on
plot(xi, X ...
        , 'Color', pen.lc ...
        , 'LineStyle', '--')

end