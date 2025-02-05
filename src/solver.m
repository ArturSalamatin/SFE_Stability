function [z,y] = solver(sigma, params, Jac, IC, A)


options = odeset(...
    'RelTol', 1e-8 ...
    , 'AbsTol', 1e-8 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
     , 'InitialStep', 1e-8 ...
    , 'MaxStep', 1e-3 ...
    , 'Jacobian', Jac ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    ..., 'Events', @(z,y) events(z,y) ...
    );

[z,y] = ode45(...
    @(z, y) my_ode(z, y, params, sigma, Jac), ...
    [A, 0], IC, options);


%% plot solutions
global DEBUG
if(DEBUG)
    for i = 1:4
        figure(700+i)
        hold on
    end
    
    mask_minus = z < 0;
    mask_plus = z >= 0;
    
    Z = z;% [exp(z(mask_minus))/2; 1-exp(z(mask_plus))/2];
    
    Y = log(abs(y));
    col = '-b';
    
    figure(701)
    plot(Z, Y(:,1), col, 'LineWidth', 1)
    %     axis([0 1 -3 1])
    
    figure(702)
    plot(Z, Y(:,2), col, 'LineWidth', 1)
    %     axis([0 1 -0.1 10])
    
    figure(703)
    plot(Z, Y(:,3), col, 'LineWidth', 1)
    %     axis([0 1 -30 10])
    
    figure(704)
    plot(Z, Y(:,4), col, 'LineWidth', 1)
    %     axis([0 1 -30 10])
end
end



