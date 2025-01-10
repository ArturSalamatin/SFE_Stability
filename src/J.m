function out = J(sigma, params)
[z,y] = solver(sigma, params);
out = log10(abs(y(end, 2))+abs(z(end)-params.z2));
end

function [z,y] = solver(sigma, params)

options = odeset(...
    'RelTol', 1e-8 ...
    , 'AbsTol', 1e-8 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
    ..., 'InitialStep', 1e-8 ...
    ..., 'MaxStep', 1e-5 ...
    , 'Jacobian',@(z, y) Jac(z, y, params, sigma) ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    , 'Events', @(z,y) events(z,y) ...
    );

eps = 1e-12;
z2 = params.z2;

[z,y] = ode15s(...
    @(z, y) my_ode(z, y, params, sigma), ...
    [0, z2-eps], [0; 0; 1], options);

z(end) - params.z2

%% plot solutions
global DEBUG
if(DEBUG)
    for i = 1:3
        figure(700+i)
        hold off
    end
    
    figure(701)
    plot(z, y(:,1), '-k', 'LineWidth', 1)
%     axis([0 1 -3 1])
    
    figure(702)
    plot(z, y(:,2), '-k', 'LineWidth', 1)
%     axis([0 1 -0.1 10])
    
    figure(703)
    plot(z, y(:,3), '-k', 'LineWidth', 1)
%     axis([0 1 -30 10])
end
end

function [value,isterminal,direction] = events(~,y)
value = [y(2); abs(y(2)) - 15];     % Detect velocity = 0
isterminal = [1; 1];   % Stop the integration
direction = [0; 0];   % Negative direction only
end



