clc
close all
clear all

global x_figid dxdt_figid c_figid
dxdt_figid = 114;
x_figid = 124;
c_figid = 134;

% T = 0.25; %linspace(1e-2,0.5, 11);
% r = 0.05;
% Sigma = -3;
% for ii = 1:numel(T)
%     t = T(ii);
%     z = linspace(0, z2(t, r));
%     x = X(t, z, r);
%     figure(x_figid)
%     hold on
%     axis([0 1 0 1])
%     plot(z/z(end), x.^2, '-k')
%     dxdt = dXdt(t, x, r);
% %     figure(dxdt_figid)
% %     hold on
% %     plot(z/z(end), dxdt)
%     
%     [z, X] = X_xi(t, r);
%     c = C(t, X, r);
%     figure(c_figid)
%     hold on
%     axis([0 1 0 1])
%     plot(z/z(end), 1-c, '-k')
%     
%     
%     figure(200)
%     hold on
%     f = 1-c + Sigma*X.^2;
%     plot(z/z(end), f)
%     
% end

t = 0.45;
r = 0.5;
R = 5;
alpha = 1*3.83;

xi_max = 1/(1-r)*log(G(sqrt(2*t), r)/ r);
sigma = (linspace((-1/(2*t)-171e-1), -1.111, 2));
val = zeros(size(sigma));
for i = 1:numel(sigma)
    i
    val(i) = J(sigma(i), t, r, R, alpha, xi_max);
    figure(400)
% hold on
    plot(sigma(1:i), val(1:i))
end

% figure(400)
% hold on
% plot(sigma, val)

% return;
% 
% t = 0.2;
% r = 0.5;

% [Sigma, val] = sigma(t, r, R, alpha);



