clc
% close all
clear all


params.t = 0.3;
params.a = sqrt(2*params.t);


params.R = 0.1;
params.alpha = 100;
params.sigma = 1;


global DEBUG
DEBUG = true;

n = 1001;
y = linspace(-20, 0, n);
eigvals = zeros(4, n);



    [V, D] = (eig(Jac_minus(-20, [], params),  eye(4), 'qz'))


for i = 1:n
    eigvals(:, i) = (eig(Jac_minus(y(i), [], params), eye(4), 'qz'));
end

figure(1)
hold on
plot(y, eigvals)



% eps = 1e-3;
% y = linspace(eps,1-eps, n);
% eigvals = zeros(4, n);
% 
% for i = 1:n
%     eigvals(:, i) = sort(eig(Jac(y(i), [], params)));
% end
% 
% figure(2)
% plot(y, eigvals)