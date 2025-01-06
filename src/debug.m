clc
clear all
close all

params.a0 = 0.02;
params.r = 0.5;
params.g1 = 1-params.r;
params.g0 = params.g1 + params.r/params.a0;

%% verify x(z,t)
figure(100)
hold on
axis([0 0.7 0 1])
box on
xlabel('z')
ylabel('X')
plot([0 0.7], params.a0*[1 1], 'k--')

t = linspace(1E-2, 0.5, 11);
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z0 = z0(params);
    params.z2 = z2(params);
    
    x = x_grid(params);
    z = z_of_x(x, params);
    
    plot(z, x, 'k-', 'LineWidth', 1)
    plot(params.z0, params.a0, 'ok', 'MarkerFaceColor', 'black')
    plot(params.z2, 0, 'sk', 'MarkerFaceColor', 'black')
end

%% verify c(z,t)
figure(200)
hold on
axis([0 0.7 0 1])
box on
xlabel('z')
ylabel('C')
% plot([0 0.7], params.a0*[1 1], 'k--')

t = linspace(1E-2, 0.5, 11);
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z0 = z0(params);
    params.z2 = z2(params);
    
    x = x_grid(params);
    z = z_of_x(x, params);
    c = c_of_x(x, params);
    
    plot(z, c, 'k-', 'LineWidth', 1)
    plot(params.z0, c_of_x(params.a0, params), 'ok', 'MarkerFaceColor', 'black')
    plot(params.z2, 1, 'sk', 'MarkerFaceColor', 'black')
end

%% verify G(x)/x
figure(300)
hold on
axis([0 2 0 params.g0])
box on
xlabel('x')
ylabel('G(x)/x')

params.t = 2;
params.a = sqrt(2*params.t);

x = x_grid(params, 0.0001);
out = G2X(x, params);

plot(x, out, 'k-', 'LineWidth', 1)
plot(params.a0*[1 1], [0 params.g0], 'k--')
plot([1 1], [0 params.g0], 'k--')


