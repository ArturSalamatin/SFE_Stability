clc
clear all
% close all


params.t = 0.3;
params.a = sqrt(2*params.t);


params.R = 0.1;
params.alpha = 100;
sigma = -0.9;


global DEBUG
DEBUG = true;
% out = J(sigma, params);


sigma = linspace(-2, -0.5, 15);
out = zeros(size(sigma));
for i = 1 :numel(sigma)
    i
    out(i) = J(sigma(i), params);
end
% 
figure(2434)
hold on
plot(sigma, out)

return
% 
% 
% [out, val] = sigma(params);
% 
% global DEBUG
% DEBUG = true;
% solver(out, params);

% figure(1000)
% plot(sigma, out)

% [out, val] = Sigma(params);



