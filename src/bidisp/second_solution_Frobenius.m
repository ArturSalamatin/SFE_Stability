function [Psi, Y] = second_solution_Frobenius(xi, alpha, tau, sigma, C1, C2, zeta2)
% coefficient derivation
%https://chat.deepseek.com/share/9nvlb0xuxegl9ixxjp

% SECOND_SOLUTION Calculate the second Frobenius solution for Psi(xi) and Y(xi)
% Inputs:
%   xi     - point where to evaluate the solution (scalar or vector)
%   alpha  - parameter
%   tau    - parameter  
%   sigma  - parameter
%   C1     - parameter
%   C2     - parameter
%   zeta2  - parameter
%
% Outputs:
%   Psi, Y - solution values at xi

%     s = sqrt(2*tau);
%     K = alpha + (1 - alpha)*s;
    r = (2 + sigma + C2) / C2;
    
    % Calculate expansion coefficients
    [A0, A1, A2, A3, B0, B1, B2, D0, D1, D2, D3, E0, E1, E2] = ...
        calculate_coefficients(alpha, tau, zeta2, sigma, C2);
    
    % Calculate series coefficients for the second solution
    [p1, p2, p3, q0, q1, q2, q3] = calculate_series_coefficients(...
        alpha, tau, sigma, C1, C2, zeta2, r, A0, A1, A2, A3, B0, B1, B2, D0, D1, D2, D3, E0, E1, E2);
    
    % Evaluate solution
    if isscalar(xi)
        [Psi, Y] = evaluate_solution(xi, r, p1, p2, p3, q0, q1, q2, q3);
    else
        Psi = zeros(size(xi));
        Y = zeros(size(xi));
        for i = 1:numel(xi)
            [Psi(i), Y(i)] = evaluate_solution(xi(i), r, p1, p2, p3, q0, q1, q2, q3);
        end
    end
end

function [A0, A1, A2, A3, B0, B1, B2, D0, D1, D2, D3, E0, E1, E2] = ...
         calculate_coefficients(alpha, tau, zeta2, sigma, C2)
    
    s = sqrt(2*tau);
    K = alpha + (1 - alpha)*s;
    
    % Leading coefficients
    A0 = 1 - alpha;
    B0 = 1 - alpha;
    D0 = 2 + sigma + C2;
    E0 = 1;
    
    % First-order coefficients
    A1 = (1 - alpha) * zeta2 * (2*alpha + (1 - alpha)*s) / s^2;
    B1 = (1 - alpha) * zeta2 * K / s^2;
    D1 = zeta2 * (2*alpha + (1 - alpha)*s) / s^2;
    E1 = zeta2 * K / s^2;
    
    % Second-order coefficients
    A2 = (1 - alpha) * zeta2^2 * (K*(3*K + alpha)/s^4 - (1 - alpha)*(alpha/2 + 2*K)/s^3);
    B2 = (1 - alpha) * zeta2^2 * (K^2 + (alpha*K)/2) / s^4;
    D2 = zeta2^2 * (K*(3*K + alpha)/s^4 - (1 - alpha)*(alpha/2 + 2*K)/s^3);
    E2 = zeta2^2 * K * (K + alpha/2) / s^4;
    
    % Third-order coefficients
    term1 = (1 - alpha)*s * (alpha*(3*alpha + 2*(1 - alpha)*s)/6 + 2*alpha*K + 3*K^2);
    term2 = K * (8*alpha^2 + (35/3)*alpha*(1 - alpha)*s + 4*(1 - alpha)^2*s^2);
    A3 = (1 - alpha) * zeta2^3 * (-term1 + term2) / s^6;
    D3 = A3 / (1 - alpha);
end

function [p1, p2, p3, q0, q1, q2, q3] = calculate_series_coefficients(...
        alpha, tau, sigma, C1, C2, zeta2, r, A0, A1, A2, A3, B0, B1, B2, D0, D1, D2, D3, E0, E1, E2)
    
    s = sqrt(2*tau);
    K = alpha + (1 - alpha)*s;
    
    % Calculate q0
    q0 = -C1 * r / (1 - alpha);
    
    % Calculate q1
    q1 = (D1 * q0 + E0) / C2;
    
    % Calculate p1
    p1 = -(A0 * q1 + A1 * q0 + B0) / (C1 * (r + 1));
    
    % Calculate q2
    q2 = (D1 * q1 + D2 * q0 + E0 * p1 + E1) / (2 * C2);
    
    % Calculate p2
    p2 = -(A0 * q2 + A1 * q1 + A2 * q0 + B0 * p1 + B1) / (C1 * (r + 2));
    
    % Calculate q3
    q3 = (D1 * q2 + D2 * q1 + D3 * q0 + E0 * p2 + E1 * p1 + E2) / (3 * C2);
    
    % Calculate p3
    p3 = -(A0 * q3 + A1 * q2 + A2 * q1 + A3 * q0 + B0 * p2 + B1 * p1 + B2) / (C1 * (r + 3));
end

function [Psi, Y] = evaluate_solution(xi, r, p1, p2, p3, q0, q1, q2, q3)
    % Evaluate the series solution up to xi^3
    if xi == 0
        Psi = 0;
        Y = 0;
    else
        xi_r = xi^r;
        Psi = xi_r * (1 + p1*xi + p2*xi^2 + p3*xi^3);
        Y = xi_r * (q0 + q1*xi + q2*xi^2 + q3*xi^3);
    end
end

% Example usage function
function example_usage()
    % Example parameters
    alpha = 0.5;
    tau = 1.0;
    sigma = 0.1;
    C1 = 1.0;
    C2 = 1.0;
    zeta2 = 0.1;
    
    % Calculate solution at multiple points
    xi_points = [0, 0.001, 0.01, 0.1, 0.5, 1.0];
    
    fprintf('Second Solution Values:\n');
    fprintf('xi\t\tPsi(xi)\t\tY(xi)\n');
    fprintf('--------------------------------\n');
    
    for i = 1:length(xi_points)
        xi = xi_points(i);
        [Psi, Y] = second_solution(xi, alpha, tau, sigma, C1, C2, zeta2);
        fprintf('%.3f\t%.6e\t%.6e\n', xi, Psi, Y);
    end
    
    % Plot the solution
    xi_plot = linspace(0, 1, 100);
    [Psi_plot, Y_plot] = second_solution(xi_plot, alpha, tau, sigma, C1, C2, zeta2);
    
    figure;
    subplot(2,1,1);
    plot(xi_plot, Psi_plot, 'b-', 'LineWidth', 2);
    xlabel('\xi');
    ylabel('\Psi(\xi)');
    title('Second Solution: \Psi(\xi)');
    grid on;
    
    subplot(2,1,2);
    plot(xi_plot, Y_plot, 'r-', 'LineWidth', 2);
    xlabel('\xi');
    ylabel('Y(\xi)');
    title('Second Solution: Y(\xi)');
    grid on;
end