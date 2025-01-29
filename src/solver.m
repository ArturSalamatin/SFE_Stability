function [z,y] = solver(sigma, params)


options = odeset(...
    'RelTol', 1e-10 ...
    , 'AbsTol', 1e-10 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
    ... , 'InitialStep', 1e-8 ...
     , 'MaxStep', 1e-3 ...
    , 'Jacobian',@(z, y) Jac(z, y, params, sigma) ...
    , 'Mass', @(z, y) Mass(z, y, params, sigma) ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    ..., 'Events', @(z,y) events(z,y) ...
    );

% a = params.a;
% z_end = a + 1/(a*sigma);
z2 = params.z2;

% if((z_end < 0) || (z_end > z2))
    z_end = z2;
% end

[z,y] = ode23tb(...
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
    plot(z/z2, y(:,1), '-k', 'LineWidth', 1)
%     axis([0 1 -3 1])
    
    figure(702)
    plot(z/z2, y(:,2), '-k', 'LineWidth', 1)
%     axis([0 1 -0.1 10])
    
    figure(703)
    plot(z/z2, y(:,3), '-k', 'LineWidth', 1)
%     axis([0 1 -30 10])
end
end

function [value,isterminal,direction] = events(~,y)
value = [y(2); abs(y(2)) - 15];     % Detect velocity = 0
isterminal = [1; 1];   % Stop the integration
direction = [0; 0];   % Negative direction only
end

function out = Mass(z, ~, params, sigma)

x = x_of_z(z, params);
dxdt = dXdt(x, params);

out = eye(3,3);
out(2,2) = (sigma*x+dxdt);

end


