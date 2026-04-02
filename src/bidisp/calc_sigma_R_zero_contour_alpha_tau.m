function [tau0, a0, out] = calc_sigma_R_zero_contour_alpha_tau(...
    solver, N, alpha, tau0, a0, pens)
global xBarLeft
R = 0;
A0 = a0;
%% calculate sigma
out = -2 + alpha./(alpha+(1-alpha).*a0);
for j = 1:size(out, 1)
    j
    guess = out(j,1);
    for i = 2:size(out, 2)
        Alpha = alpha(j,i);
        Tau0 = tau0(j,i);
        
        params = poly_case(A0, Alpha, Tau0, R);
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
