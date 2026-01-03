function [tau0, a0, out] = calc_sigma_full(...
    solver, N, alpha, tau0, a0, R_vals, h_vals, pens)
%% calculate sigma
out = -2 + alpha./(alpha+(1-alpha).*a0) ...
    + zeros(numel(R_vals), numel(h_vals));
for i = 1:numel(h_vals)
    disp(['i = ', num2str(i)]);
    h = h_vals(i);
    guess = out(1,i);
    for j = 1:numel(R_vals)
        R = R_vals(j);
        params = poly_case(a0, alpha, tau0, R, h);
        params.marker = 's';
        params.pen = set_pen(...
            pens.lc{min(1, numel(pens.lc))}, ...
            pens.style{min(j, numel(pens.style))});
        mesh = set_full_mesh(N, params, 0);
        starter = @(sigma, params) starter_full(...
            @(problem, mesh)solver(problem, mesh, params), sigma, params, mesh);
        % dust fraction has been extracted at the inlet
        sigma = fit_sigma(starter, params, guess);
        %% save sigma to the output matrix
        out(j,i) = sigma;
        %% remeber last sigma as a guess
        guess = sigma;
    end
end