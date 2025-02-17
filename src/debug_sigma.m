
fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)

clc
clear all
% close all

global A fig_id
A = -35;
fig_id = 700;

params.t = 0.4;
params.a = sqrt(2*params.t);

params.R = 0.2;
params.alpha = 1;

% sigma = 0.5;
% % global DEBUG
% % DEBUG = true;


for i = [1:5]
    figure(fig_id+i)
    legend
    box on
    hold off
end

Sigma = linspace(-15, 15, 15);
Sigma = 0.1;
f = zeros(size(Sigma));
for i = 1:numel(Sigma)
    left_solver(params, Sigma(i));
    %     right_solver(params, Sigma(i));
end

figure(fig_id+1)
xlabel(['{\xi}'])
ylabel(['{\Phi}'])

figure(fig_id+2)
xlabel(['{\xi}'])
ylabel(['{\Psi}'])

figure(fig_id+3)
xlabel(['{\xi}'])
ylabel(['X'])

figure(fig_id+4)
xlabel(['{\xi}'])
ylabel(['{\Gamma}'])
for i = [1:5]
    figure(fig_id+i)
    legend
end

figure(8034)
hold on
plot(Sigma, f)

% A = -15;
% for i = 1:numel(Sigma)
%     right_solver(params, Sigma(i));
% end

% A = -35;
% out = J(0.1, params);
% A = -30;
% out = J(0.1, params);
% A = -25;
% out = J(0.1, params);
% A = -20;
% out = J(0.1, params);
% A = -15;
% out = J(0.1, params);
% A = -10;
% out = J(0.1, params);


%
% sigma = linspace(0.001, 2.5, 45);
% out = zeros(size(sigma));
% for i = 1 :numel(sigma)
%     i
%     out(i) = J([sigma(i), 0, 0, -1], params);
% end
% %
% figure(2434)
% hold on
% plot(sigma, out)
%
return


[out, val] = sigma(params);

global DEBUG
DEBUG = true;
% J([0.320940754176750,-0.372695959119099,3.816116788745306], params);
% J([0.334339751816437,-0.371433948805199,3.817722879358429], params);
% J([0.135703653179201,-0.391489199199079,3.791934999365223], params);

J(out, params);

% figure(1000)
% plot(sigma, out)

% [out, val] = Sigma(params);
%
% global DEBUG
% DEBUG = true;
% J(out, params);


