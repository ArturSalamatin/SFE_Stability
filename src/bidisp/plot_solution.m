function plot_solution(...
    t, y, pen, params, sol)
global xBarLeft
fig_id = 700;

factor = 1;

names = {'{\Psi}','X','{\Phi}','{\Gamma}','{\Omega}','Y', '{\Psi}+X', 'Q', 'P', 'P^{\prime}'};
plot_solution_assymptotics(fig_id, pen, params, sol);
x_left = linspace(1,1-2*xBarLeft,10001)*params.a;
[~,xi, Psi, X] = calc_inlet_solution_assymptotics(...
x_left, params, sol.sigma ...
...fig_id,pen, params, sol
    );
Psi = Psi/sol.factor;
X = X/sol.factor;

figure(fig_id+1)
hold on
plot(xi, Psi ...
        , 'Color', 'k' ...
        , 'LineStyle', '-.')

figure(fig_id+2)
hold on
plot(xi, X ...
        , 'Color', 'k' ...
        , 'LineStyle', '-.')

my_figure(fig_id+10)
hold on
    plot(xi, Psi./X ...
        , 'Color', 'k' ...
        , 'LineStyle', '-.')




a0 = params.a0;
% x = a0;% linspace(a0, 2*a0, 101);
% xi = z_of_x(x, params)/sol.base_state.z2;
xi0 = z_of_x(a0, params)/params.z2;
g1 = params.g1;
g0 = params.g0 - g1;
a = params.a;
dz0dt = params.dz0dt;
C1 = params.C1;
C2 = params.C2;
sigma = sol.sigma;
psi0 = g0*a*dz0dt - (g0+g1)/C1*(C2 + (1+sigma)*(1-xi0));

figure(fig_id+1)
hold on
plot(xi0, psi0, 'ks', 'markerfacecolor', 'k')

for i = [1,2] %,5,6]% 1:6%numel(names)
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
    if(i == 2)
        axis([0 1 0 1])
    end
end

my_figure(fig_id+10)
hold on
plot(t, y(:,1)./y(:,2)...
        , 'LineWidth', 1 ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-'... pen.style{1}(1) ...
        )


% figure(fig_id+8)
% axis([-Inf Inf -Inf 5])

figure(701)
end