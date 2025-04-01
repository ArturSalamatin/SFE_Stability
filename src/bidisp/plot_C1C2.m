clc
close all

params.a0 = 0.1;
params.r = 0.1;
params.a1 = 1;
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

t = linspace(0, 0.5, 301);
I1 = round(numel(t)/4.5);
I2 = round(3*numel(t)/4);
params.t = t;
params.a = sqrt(2*params.t);
a = params.a;



h1 = my_figure(300);
hold on
xlabel('\tau')
ylabel('{\itC}_1')

h2 = my_figure(400);
hold on
xlabel('\tau')
ylabel('{\itC}_2')

alpha = [0.1, 0.3, 0.5, 0.7, 1];
alpha = [0.1, 0.3, 0.5, 0.7];
for i = 1:numel(alpha)
    
params.r = alpha(i);
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

C1 = a./z2(params);
C2 = C1.*a.*dz2dt(params);

figure(300)
plot(t, C1, 'k-', 'LineWidth', 1)
plot(t(I1), C1(I1), 'ok', 'MarkerFaceColor', 'black')

figure(400)
plot(t, C2, 'k-', 'LineWidth', 1)
plot(t(I1), C2(I1), 'ok', 'MarkerFaceColor', 'black')
end



params.r = 0.5;
a0 = [0.1, 0.3, 0.5];
for i = 1:numel(a0)
    params.a0 = a0(i);
    
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

C1 = a./z2(params);
C2 = C1.*a.*dz2dt(params);

figure(300)
plot(t, C1, 'k--', 'LineWidth', 1)
plot(t(I2), C1(I2), 'sk', 'MarkerFaceColor', 'black')

figure(400)
plot(t, C2, 'k--', 'LineWidth', 1)
plot(t(I2), C2(I2), 'sk', 'MarkerFaceColor', 'black')
end

saveas(h1, 'Figs\C1.emf')
saveas(h1, 'Figs\C1.fig')
saveas(h2, 'Figs\C2.emf')
saveas(h2, 'Figs\C2.fig')

