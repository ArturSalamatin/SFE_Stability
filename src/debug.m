clc
clear all
close all

params.a0 = 0.05;
params.r = 0.5;
params.g1 = 1-params.r;
params.g0 = params.g1 + params.r/params.a0;


figure(100)
hold on
axis([0 0.7 0 1])
box on
xlabel('z')
ylabel('X')

t = linspace(1E-2, 0.5, 11);
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z0 = z0(params);
    params.z2 = z2(params);
    
    x = x_grid(params);
    z = z_of_x(x, params);
    
    plot(z, x, 'k-', 'LineWidth', 1)
end


plot([0 0.7], params.a0*[1 1], 'k--')



