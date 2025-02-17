function out = J(x, params)
global A
out = 0;

R = params.R;
alpha = params.alpha;
alpha2 = (alpha)^2;
a = params.a;
a2 = a*a;

sigma = x(1);
%%
[t, val] = left_solver(params, sigma);

%%
phi0 = x(2);
psi0 = x(3);
chi0 = -psi0;
gamma0 = -a*alpha*phi0;

c0 = [phi0, psi0, chi0, gamma0];

h1 = exp(-abs(A))/2;
phi1 = -gamma0;
chi1 = -phi0;
psi1 = 2*phi0-(1+sigma)*chi0;
gamma1 = R*gamma0 - a2*alpha2*(phi0+R*psi0);
c1 = [phi1, psi1, chi1, gamma1];

h2 = h1*h1;
phi2 = -gamma1/2;
chi2 = (gamma0 + (1+sigma)*(phi0-chi0))/4;
psi2 = chi2-gamma0;
gamma2 = (R*gamma1 - a2*alpha2*(phi1+R*psi1))/2;
c2 = [phi2, psi2, chi2, gamma2];

IC = c0 + c1*h1 + c2*h2;
    
Jac = @(z, u) Jac_plus(z, u, params, sigma);
[z_plus, u_plus] = solver(sigma, params, Jac, IC, [abs(A), 0]);
plot_solution(1-exp(-z_plus)/2,u_plus,'-k');
%%
% c1 = exp(-abs(A))/2;
% phi1 = 1;
% psi1 = 0;
% chi1 = 0;
% gamma1 = -R;
% 
% c2 = c1*c1;
% phi2 = gamma1/2;
% psi2 = -phi1/2; 
% chi2 = psi2/(-sigma)+chi1; % phi1/(2*sigma);
% gamma2 = -(R*gamma1 - a2*alpha2*(phi1+R*psi1))/2;
% IC = ([0;0;0;1] + ...
%     [phi1; psi1; chi1; gamma1]*c1 + ...
%     [phi2; psi2; chi2; gamma2]*c2)';
% Jac = @(z, u) Jac_minus(z, u, params, sigma);
% [z_minus, u_minus] = ...
%     solver(sigma, params, Jac, IC, [-abs(A), 0]);

% c1 = exp(-AA)/2;
% phi1 = 1;
% psi1 = 0;
% chi1 = 0;
% gamma1 = -R;
% 
% c2 = c1*c1;
% phi2 = gamma1/2;
% psi2 = -1/2; %-(phi1+chi1)/2;
% chi2 = - psi2/sigma; % phi1/(2*sigma);
% gamma2 = -(R*gamma1 - a2*alpha2*(phi1+R*psi1))/2;
% 
% BC = ([0;0;0;1] + ...
%     [phi1; psi1; chi1; gamma1]*c1 + ...
%     [phi2; psi2; chi2; gamma2]*c2)';
% 
% FACTOR = u_minus(end, 4)/BC(4);
% BC = FACTOR*BC;

%%
out = log10...
    (sum(abs(u_plus(end, :)-val)));

out
x

% 
% z = [exp(z_minus)/2; 1-exp(-z_plus(end:-1:1))/2];
% u = [u_minus;        u_plus(end:-1:1, :)];

end

function plot_solution(t,y,col)
global fig_id
for i = 1:4
    figure(fig_id+i)
    axis([0 1 -Inf Inf])
    plot(t,y(:,i),col, 'LineWidth', 1)
    hold on
end
end




