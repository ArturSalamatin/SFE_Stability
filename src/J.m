function out = J(sigma, params)
[z,y] = solver(sigma, params);
out = log10(abs(y(end, 2))+abs(z(end)-params.z2));
end

function [z,y] = solver(sigma, params)

options = odeset(...
    'RelTol', 1e-6 ...
    , 'AbsTol', 1e-6 ...
    , 'NormControl', 'on' ...
    ..., 'NonNegative', [1,3] ...
    ..., 'InitialStep', 1e-8 ...
    ..., 'MaxStep', 1e-5 ...
    , 'Jacobian',@(z, y) Jac(z, y, params, sigma) ...
    ..., 'Stats','on' ...
    ...,'OutputFcn', @odeplot ...
    );

z0 = params.z0;
% takes into account jump at z = z0
eps = 1e-12;

[z,y] = ode45(...
    @(z, y) my_ode(z, y, params, sigma), ...
    [0, z0-eps], [0; 0; 1], options);

a0 = params.a0;
r = params.r;
dZ0dt = dz0dt(params);
factor = (2*a0*a0+r*dZ0dt)/(2*a0*a0-r*dZ0dt);
z2 = params.z2;

options = odeset(options, 'Events', @(z,y) events(z,y));
[Z,Y,te,ye,ie] = ode45(...
    @(z, y) my_ode(z, y, params, sigma), ...
    [z0+eps, z2-1e-6], y(end, :).*[1,factor,1], options);
Z(end)

%% plot solutions
global DEBUG
if(DEBUG)
    figure(701)
    hold off
    plot(z, y(:,1), '-k', 'LineWidth', 1)
    hold on
    plot(Z, Y(:,1), '-r', 'LineWidth', 1)
%     axis([0 1 -3 1])
    
    figure(702)
    hold off
    plot(z, y(:,2), '-k', 'LineWidth', 1)
    hold on
    plot(Z, Y(:,2), '-r', 'LineWidth', 1)
%     axis([0 1 -0.1 10])
    
    figure(703)
    hold off
    plot(z, y(:,3), '-k', 'LineWidth', 1)
    hold on
    plot(Z, Y(:,3), '-r', 'LineWidth', 1)
%     axis([0 1 -30 10])
end

z = [z;Z];
y = [y;Y];
end


function out = chi(z, psi, sigma, params)
X0 = x_of_z(z, params);
dx0dt = dXdt(X0, params);

out = -psi./(sigma*X0 + dx0dt);

end

function [value,isterminal,direction] = events(~,y)
value = y(2);     % Detect velocity = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end



