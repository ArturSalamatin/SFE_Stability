clc
clear all
close all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)

params.a0 = 0.3;
params.a1 = 1.0;
params.r = 0.25; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

t = [linspace(1E-8, 0.5, 6), (params.a0^2)/2];
t = sort(t);

%% verify x(z,t)
verify_x(params, t);
%% verify c(z,t)
verify_c(params, t);

%% verify G(x)/x
% verify_GdivX(params);

saveas(200, 'Figs/Fig2C-R1.eps');
saveas(200, 'Figs/Fig2C-R1.emf');
saveas(100, 'Figs/Fig2D-R1.eps');
saveas(100, 'Figs/Fig2D-R1.emf');
saveas(200, 'Figs/Fig2C-R1.fig');
saveas(100, 'Figs/Fig2D-R1.fig');




function verify_x(params, t)
my_figure(100)
hold on
axis([0 0.8 0 1])
% box on
xlabel(['\zeta, ' char(8211)])
ylabel(['{x}_0, ' char(8211)])
box on

set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])



plot([0 1], params.a0*[1 1], 'k--')


my_figure(110)
hold on
axis([0 1 0 1])
% box on
xlabel('\xi')
ylabel('$\bar{x}_0$', 'interpreter', 'latex')
% plot([0 0.7], params.a0*[1 1], 'k--')

for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z0 = z0(params);
    params.z0
    params.z2 = z2(params);
    
    x = x_grid(params);
    z = z_of_x(x, params);
    
    params.x = x;
    params.z = z;
    
    %     x_inv = x_of_z(z, params);
    
    figure(100)
    if(i==2)
    plot(z, x, 'b-', 'LineWidth', 1)
    else
    plot(z, x, 'k-', 'LineWidth', 1)
    end
    %     plot(z, x_inv, 'r--', 'LineWidth', 2)
    if(params.z0>=0)
        plot(params.z0, params.a0, 'ok', 'MarkerFaceColor', 'black')
    end
    plot(params.z2, 0, 'sk', 'MarkerFaceColor', 'black')
    
    figure(110)
    if(i==2)
    plot(z/params.z2, x/params.a, 'b-', 'LineWidth', 1)
    else
    plot(z/params.z2, x/params.a, 'k-', 'LineWidth', 1)
    end
    %     plot(z, x_inv, 'r--', 'LineWidth', 2)
    plot(params.z0/params.z2, params.a0/params.a, 'ok', 'MarkerFaceColor', 'black')
    %     plot(params.z2, 0, 'sk', 'MarkerFaceColor', 'black')
end
figure(100)
annotation('textbox',...
    [0.783214285714285 0.783285712741676 0.082142858675548 0.0885714301154727],...
    'String',{'D'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
end

function verify_c(params, t)
my_figure(200)
hold on
axis([0 0.8 0 1])
xlabel(['\zeta, ' char(8211)])
ylabel(['{c}_0, ' char(8211)])

box on

set(gcf, 'units', 'centimeters', 'OuterPosition', [10.42 6.27 8.5 9])
set(gca, 'FontSize', 10, 'Position', [0.16, 0.2, 0.78, 0.74])


my_figure(210)
hold on
axis([0 1 0 1])
xlabel(['\zeta, ' char(8211)])
ylabel(['{c}_0, ' char(8211)])

% t = linspace(1E-2, 0.5, 11);
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z0 = z0(params);
    params.z2 = z2(params);
    
    x = x_grid(params);
    z = z_of_x(x, params);
    c = c_of_x(x, params);
    
    figure(200)
    if(i==2)
    plot(z, c, 'b-', 'LineWidth', 1)
    else
    plot(z, c, 'k-', 'LineWidth', 1)
    end
    plot(params.z0, c_of_x(params.a0, params), 'ok', 'MarkerFaceColor', 'black')
    plot(params.z2, 1, 'sk', 'MarkerFaceColor', 'black')
    
    figure(210)
    if(i==2)
    plot(z/params.z2, c, 'b-', 'LineWidth', 1)
    else
    plot(z/params.z2, c, 'k-', 'LineWidth', 1)
    end
    plot(params.z0/params.z2, c_of_x(params.a0, params), 'ok', 'MarkerFaceColor', 'black')
    %     plot(params.z2, 1, 'sk', 'MarkerFaceColor', 'black')
end


t = linspace((params.a0^2)/2, 0.5, 101);
z = zeros(size(t));
c = zeros(size(t));
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*t(i));
    z(i) = z0(params);
    c(i) = c_of_x(params.a0, params);
end

figure(200)
plot(z,c, 'k--', 'LineWidth', 1)
annotation('textbox',...
    [0.783214285714285 0.185190474646433 0.082142858675548 0.0885714301154727],...
    'String',{'C'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName','Times New Roman');
end

function verify_GdivX(params)
my_figure(300)
hold on
axis([0 2 0 params.g0])
xlabel('{\itx}')
ylabel('{\itG}({\itx})/{\itx}')


a0 = [0.1, 0.3, 0.5, 0.7, 0.9];

for i = 1:numel(a0)
    params.a0 = a0(i);
    params.g0 = params.g1 + params.r/params.a0;
    params.t = 2;
    params.a = sqrt(2*params.t);
    
    x = x_grid(params, 0.0001);
    out = G_div_X(x, params);
    
    plot(x, out, 'k-', 'LineWidth', 1)
    plot(params.a0*[1 1], [0 params.g0], 'k--')
    plot([1 1]*params.a1, [0 params.g0], 'k--')
end
end
