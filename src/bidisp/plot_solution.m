function plot_solution(...
    pen, params, sol)
global xBarLeft
fig_id = 700;

t = sol.t;
y = sol.y;

names = {'{\Psi}','X','{\Phi}','{\Gamma}','{\Omega}','Y', '{\Psi}+X', 'Q', 'P', 'P^{\prime}'};
plot_solution_assymptotics(fig_id, pen, params, sol.sigma);
x_left = linspace(1,1-xBarLeft*2,1001)*params.a;
[~,xi, Psi, X] = calc_inlet_solution_assymptotics(...
    x_left, params, sol.sigma );

Psi = Psi/sol.factor;
X = X/sol.factor;

figure(fig_id+1)
hold on
plot(xi, Psi ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-.')

figure(fig_id+2)
hold on
plot(xi, X ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-.')

my_figure(fig_id+10)
hold on
    plot(xi, Psi./X ...
        , 'Color', pen.lc ...
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
xi0 = params.xi0;
psi0 = -(g1+g0*(1-xi0))*a*dz0dt - (1+sigma)*a0/a;

figure(fig_id+1)
hold on
plot(xi0, psi0, 'ks', 'markerfacecolor', 'k')

factor = 1;
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

C1 = params.C1;
C2 = params.C2;
alpha = params.r;
factor = (1-alpha)*C2/(C1*(2+sigma+C2));

my_figure(fig_id+10)
hold on
plot(t, y(:,1)./y(:,2)./(t')/factor...
        , 'LineWidth', 1 ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-'... pen.style{1}(1) ...
        )


% figure(fig_id+8)
% axis([-Inf Inf -Inf 5])

my_figure(fig_id+20)
hold on
plot(y(:,1), y(:,2)...
        , 'LineWidth', 1 ...
        , 'Color', pen.lc ...
        , 'LineStyle', '-'... pen.style{1}(1) ...
        )
    xlabel(names{1})
    ylabel(names{2})

figure(701)
end