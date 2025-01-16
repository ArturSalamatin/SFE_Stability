clc
clear all
% close all


params.t = 0.45;
params.a = sqrt(2*params.t);
params.z2 = z2(params);


params.R = 1;
params.alpha = 5;
% sigma = 100;


% sigma = linspace(-30, -20, 10);
% psi0 = 5;
% out = zeros(size(sigma));
% for i = 1 :numel(sigma)
%     out(i) = J([sigma(i), psi0], params);
% end

[out, val] = sigma(params);

global DEBUG
DEBUG = true;
solver(out, params);

% figure(1000)
% plot(sigma, out)