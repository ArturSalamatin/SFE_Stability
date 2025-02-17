function [z,y] = solver(sigma, params, Jac, IC, segment)


options = odeset(...
    'RelTol', 1e-10 ...
    , 'AbsTol', 1e-10 ...
    , 'NormControl', 'on' ...
     , 'InitialStep', 1e-8 ...
    ... , 'MaxStep', 1e-3 ...
    , 'Jacobian', Jac ...
    ..., 'Stats','on' ...
    ... ,'OutputFcn', @odeplot ...
    ..., 'Events', @(z,y) events(z,y) ...
    );

[z,y] = ode15s(...
    @(z, y) my_ode(z, y, params, sigma, Jac), ...
    segment, IC, options);


%% plot solutions
% global DEBUG
% if(DEBUG)
% %     for i = 1:4
% %         figure(700+i)
% %         hold on
% %     end
%     
%     mask_minus = z < 0;
%     mask_plus = z >= 0;
%     
%     Z = [1-exp(-z(mask_plus))/2; exp(z(mask_minus))/2];
%     
%     Y = y;
%     col = '-m';
%     
%     figure(701)
%     plot(Z, Y(:,1), col, 'LineWidth', 1)
%     hold on
%     %     axis([0 1 -3 1])
%     
%     figure(702)
%     plot(Z, Y(:,2), col, 'LineWidth', 1)
%     hold on
%     %     axis([0 1 -0.1 10])
%     
%     figure(703)
%     plot(Z, Y(:,3), col, 'LineWidth', 1)
%     hold on
%     %     axis([0 1 -30 10])
%     
%     figure(704)
%     plot(Z, Y(:,4), col, 'LineWidth', 1)
%     hold on
%     %     axis([0 1 -30 10])
% end
end



