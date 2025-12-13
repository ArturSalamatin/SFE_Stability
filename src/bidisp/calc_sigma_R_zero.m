function [tau0, a0, out] = calc_sigma_R_zero(...
    solver, N, alpha, tau0, a0, pens)
global xBarLeft
R = 0;
%% calculate sigma
out = -2 + alpha./(alpha+(1-alpha).*a0);
guess = out(1,1);
for j = 1:size(out, 1)
    j
    for i = 1:size(out, 2)
        A0 = a0(j,i);
        Tau0 = tau0(j,i);
        
        if(A0*A0/2 >= Tau0)
            % dust fraction is not extrcted yet
            out(j,i) = NaN;
            guess = -1 + (1-alpha*A0)/(alpha+(1-alpha)*A0);
        else
            params = poly_case(A0, alpha, Tau0, R);
            params.marker = 's';
            params.pen = set_pen(...
                pens.lc{min(1, numel(pens.lc))}, ...
                pens.style{min(j, numel(pens.style))});
            mesh = set_left_mesh(N, params, xBarLeft);
            starter = @(sigma, params) starter_R_zero_X_Psi(...
                @(problem, mesh)solver(problem, mesh, params), sigma, params, mesh);
            % dust fraction has been extracted at the inlet
            sigma = fit_sigma(starter, params, guess);
            %% save sigma to the output matrix
            out(j,i) = sigma;
            %% remeber last sigma as a guess
            guess = sigma;
        end
    end
end
%% correction for the moment tau = (a0^2)/2
for j = 1:size(out, 1)
    for i = size(out, 2):-1:1
        A0 = a0(j,i);
        Tau0 = tau0(j,i);
        if(A0*A0/2 > Tau0)
            a0(j,i) = sqrt(2*Tau0);
            out(j,i) = -2 + alpha/(alpha+(1-alpha)*a0(j,i));
            break;
        end        
    end    
end