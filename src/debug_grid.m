clc
clear all
% close all

global A
A = -35;

params.t = 0.1;
params.a = sqrt(2*params.t);


params.R = 0.2;
params.alpha = 50;



n = 101;
m = 101;
Sigma = linspace(-2,2,n);
X2 = linspace(-10, 10, m)-10;

[Sigma, X2] = meshgrid(Sigma, X2);

out = zeros(size(X2));

for i = 1:n
    for j = 1:m
        out(j,i) = J([Sigma(j,i), X2(j,i)], params);
    end
end

figure(1)
hold on
mesh(Sigma, X2, out)





% figure(1000)
% plot(sigma, out)

% [out, val] = Sigma(params);
% 
% global DEBUG
% DEBUG = true;
% J(out, params);


