function [z,y] = solver(xx, params, h)

if(nargin == 2)
    h = 1E-8;
end

sigma = xx(1);
options = odeset(...
    'RelTol', 1e-8 ...
    , 'AbsTol', 1e-10 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
    ... , 'InitialStep', 1e-8 ...
    ... , 'MaxStep', 1e-5 ...
    , 'Jacobian',@(z, y) Jac(z, y, params, sigma) ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    ..., 'Events', @(z,y) events(z,y) ...
    );

z_end = params.z2 - h;

[z,y] = ode15s(...
    @(z, y) my_ode(z, y, params, sigma), ...
    [0, z_end], [0;0;1], options);

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

