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
params.r = 0.05; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

t = linspace(1E-2, 0.5, 5);

%% verify x(z,t)
verify_x(params, t);
%% verify c(z,t)
verify_c(params, t);
%% verify G(x)/x
% verify_GdivX(params);

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
    plot(z, x, 'k-', 'LineWidth', 1)
    %     plot(z, x_inv, 'r--', 'LineWidth', 2)
    if(params.z0>0)
        plot(params.z0, params.a0, 'ok', 'MarkerFaceColor', 'black')
    end
    plot(params.z2, 0, 'sk', 'MarkerFaceColor', 'black')
    
    figure(110)
    plot(z/params.z2, x/params.a, 'k-', 'LineWidth', 1)
    %     plot(z, x_inv, 'r--', 'LineWidth', 2)
    plot(params.z0/params.z2, params.a0/params.a, 'ok', 'MarkerFaceColor', 'black')
    %     plot(params.z2, 0, 'sk', 'MarkerFaceColor', 'black')
end
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
    plot(z, c, 'k-', 'LineWidth', 1)
    plot(params.z0, c_of_x(params.a0, params), 'ok', 'MarkerFaceColor', 'black')
    plot(params.z2, 1, 'sk', 'MarkerFaceColor', 'black')
    
    figure(210)
    plot(z/params.z2, c, 'k-', 'LineWidth', 1)
    plot(params.z0/params.z2, c_of_x(params.a0, params), 'ok', 'MarkerFaceColor', 'black')
    %     plot(params.z2, 1, 'sk', 'MarkerFaceColor', 'black')
end
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
