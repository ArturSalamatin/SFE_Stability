clc
close all

params.a0 = 0.1;
params.r = 0.1;
params.a1 = 1;
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

t = linspace(0, 0.5, 3001);
% I1 = round(numel(t)/4.5);
% I2 = round(3*numel(t)/4);
params.t = t;
params.a = sqrt(2*params.t);
a = params.a;



h5 = my_figure(500);
hold on
axis([0 0.5 0 Inf])
xlabel('\tau')
ylabel('{\zeta}_2')

alpha = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9];
for i = 1:numel(alpha)
    
params.r = alpha(i);
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

zeta2 = z2(params);

figure(500)
plot(t, zeta2, 'k-', 'LineWidth', 1)
% plot(t(I1), C1(I1), 'ok', 'MarkerFaceColor', 'black')
end



params.r = 0.5;
a0 = [0.1, 0.3, 0.5];
for i = 1:numel(a0)
    params.a0 = a0(i);
    
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

zeta2 = z2(params);

figure(500)
% plot(t, zeta2, 'k--', 'LineWidth', 1)
% plot(t(I2), C1(I2), 'sk', 'MarkerFaceColor', 'black')
end

saveas(h5, 'Figs\z2.emf')
saveas(h5, 'Figs\z2.fig')

