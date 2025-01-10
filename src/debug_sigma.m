clc
clear all
% close all

global DEBUG
DEBUG = true;

params.t = 0.45;
params.a = sqrt(2*params.t);
params.z2 = z2(params);


params.R = 1;
params.alpha = 5;
% sigma = 100;


sigma = linspace(-40, -0.001, 205);
out = zeros(size(sigma));
for i = 1:numel(sigma)
    out(i) = J(sigma(i), params);
end

figure(1000)
plot(sigma, out)