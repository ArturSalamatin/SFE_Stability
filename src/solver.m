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

z_end = 1;

[z,y] = ode23t(...
    @(z, y) my_ode(z, y, params, sigma), ...
    [1e-8, z_end], [0;1e-10;0;1], options);

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
    
    figure(704)
    plot(z, y(:,4), '-k', 'LineWidth', 1)
    %     axis([0 1 -30 10])
end
end

function out = Mass(z, ~, params, sigma)

out = eye(4,4);
out(2,3) = z;
out(3,3) = z*(1-z);
end


