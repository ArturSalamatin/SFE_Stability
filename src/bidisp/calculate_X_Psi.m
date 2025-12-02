function [Psi, X] = calculate_X_Psi(...
    xi, tau, alpha, C1, C2, Z2, sigma)
% Calculate X and Psi at a given xi > 0
% Inputs:
%   xi    - scalar value > 0
%   tau, alpha, C1, C2, Z2, sigma - parameters
%   b0    - arbitrary constant from Frobenius expansion
% Outputs:
%   X, Psi, x0, c0 - values at the given xi

b0 = -1;
% 1. First find x0 corresponding to this xi
%     x0_guess = sqrt(2*tau); % Initial guess
%     options = optimset('Display', 'off', 'TolX', 1e-12);

% Solve for x0 using equation (3)
%     x0 = fzero(@(x0) xi_equation(x0, xi, tau, alpha, Z2), x0_guess, options);

% 2. Calculate c0 from equation (4)
%     D = alpha + (1-alpha)*sqrt(2*tau);
%     c0 = 1 - (alpha + (1-alpha)*x0)/D;

% 3. Set up and solve the ODE system

% Integrate from small epsilon to xi
epsilon = xi;

% Initial conditions from Frobenius expansion at epsilon
t = (2 + sigma)/C2;
A10 = (1 - alpha)/C1;
a = sqrt(2*tau);

b1 = b0*(...
    -Z2/(2*tau)*(2*alpha+(1-alpha)*a)...
    +A10/(t+1));

% d1 = Z2*alpha/(2*tau);
% x1 = -Z2*(alpha+(1-alpha)*a)/a;
% 
% b1 = b0/C2*(d1+x1*(1+sigma)/a...
%      -(1-alpha)/(C1*(t+1)));

X = b0 * epsilon^t + b1*epsilon^(t + 1);
Psi = -A10/(t + 1) * b0 * epsilon^(t + 1);
end

function F = xi_equation(x0, xi_target, tau, alpha, Z2)
% Equation (3) for xi as function of x0
D = alpha + (1-alpha)*sqrt(2*tau);
F = (1/((1-alpha)*Z2)) * (sqrt(2*tau) - x0 + ...
    (alpha/(1-alpha)) * log((alpha + (1-alpha)*x0)/D)) - xi_target;
end

function dYdxi = ode_system(xi, Y, tau, alpha, C1, C2, Z2, sigma)
% ODE system for Psi and X
% Y = [Psi; X]

Psi = Y(1);
X = Y(2);

% Find x0 at current xi
x0_guess = sqrt(2*tau);
options = optimset('Display', 'none', 'TolX', 1e-12);
x0 = fzero(@(x0) xi_equation(x0, xi, tau, alpha, Z2), x0_guess, options);

% Calculate coefficients
D = alpha + (1-alpha)*sqrt(2*tau);
A1 = (2*tau*(1-alpha)/(C1*x0^2*D)) * (alpha + (1-alpha)*x0);
A2 = sqrt(2*tau)*(1-alpha)/(C1*x0);
B1 = (2*tau/D) * (alpha/x0 + (1-alpha)) + x0*(1 + sigma);

% ODEs
dPsidxi = -A1 * X - A2 * Psi;
dXdxi = (B1 * X + sqrt(2*tau) * Psi) / (x0 * C2 * xi);

dYdxi = [dPsidxi; dXdxi];
end

% Example usage function
function example_usage()
% Parameters (example values - adjust as needed)
tau = 1.0;
alpha = 0.3;
C1 = 1.0;
C2 = 1.0;
Z2 = 1.0;
sigma = 0.5;
b0 = 1.0; % arbitrary constant

% Calculate at specific xi values
xi_values = [0.1, 0.5, 1.0, 2.0];

fprintf('X and Psi values for different xi:\n');
fprintf('xi\t\tX\t\tPsi\t\tx0\t\tc0\n');

for i = 1:length(xi_values)
    xi = xi_values(i);
    [X, Psi, x0, c0] = calculate_X_Psi(xi, tau, alpha, C1, C2, Z2, sigma, b0);
    fprintf('%.1f\t%.6f\t%.6f\t%.6f\t%.6f\n', xi, X, Psi, x0, c0);
end
end