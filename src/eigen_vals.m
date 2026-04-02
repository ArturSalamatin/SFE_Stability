
fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)

clc
close all
clear all


params.t = 0.3;
params.a = sqrt(2*params.t);


params.R = 0.3;
params.alpha = 30;
params.sigma = 0.2;


global DEBUG
DEBUG = true;
%%
n = 100701;
eps = 1e-7;
y = linspace(eps, 1-eps, n);
eigvals = zeros(4, n);
sigma = params.sigma;
for i = 1:n
    eigvals(:, i) = ...
    (eig(Jac_changed_vars(y(i), [], params, sigma)));
end

figure(1)
hold on
box on
axis([0 1 -250 250])
h1 = ...
    plot(y, imag(eigvals(1,:)), 'r.', 'DisplayName', 'imag');
h2 = ...
    plot(y, real(eigvals(1,:)), 'k.', 'DisplayName', 'real');

plot(y, imag(eigvals), 'r.')
plot(y, real(eigvals), 'k.')
legend([h1, h2])
xlabel(['{\xi}'])
ylabel(['{eigenvals}'])
title(['r = 0.3, f = 800, \sigma = 0.2'])
return

%%
n = 10001;
y = linspace(-15, 15, n);
eigvals = zeros(4, n);
sigma = params.sigma;
for i = 1:n
    eigvals(:, i) = ...
    (eig(Jac_changed_vars(y(i), [], sigma, params)));
end

figure(1)
hold on
% axis([-Inf Inf -250 250])
plot(y, real(eigvals), 'k.')
plot(y, imag(eigvals), 'r.')

return

%%
n = 10001;
y = linspace(-5, 1, n);
eigvals = zeros(4, n);

for i = 1:n
    eigvals(:, i) = ...
    (eig(Jac_minus(y(i), [], params)));
end

figure(1)
hold on
axis([-Inf Inf -250 250])
plot(y, real(eigvals), 'k.')
% plot(y, imag(eigvals), 'r.')




%%
y = linspace(5, -1, n);
for i = 1:n
    eigvals(:, i) = ...
    (eig(Jac_plus(y(i), [], params)));
end

figure(1)
hold on
% axis([-Inf Inf -50 50])
plot(y, real(eigvals), 'b.')
% plot(y, imag(eigvals), 'g.')


%%
% eps = 1e-4;
% y = linspace(eps,1-eps, n);
% eigvals = zeros(4, n);
% 
% for i = 1:n
%     eigvals(:, i) = eig(Jac(y(i), [], params));
% end
% 
% figure(2)
% hold on
% axis([0 1 -100 100])
% plot(y, real(eigvals), 'k.')
% % plot(y, imag(eigvals), 'r.')