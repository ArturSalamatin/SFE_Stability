function out = calc_sigma_full_2(...
    solver, N, alpha, Tau0_vals, a0, R, h_vals, pens)
%% calculate sigma
out = -2 + alpha./(alpha+(1-alpha).*a0) ...
    + zeros(numel(Tau0_vals), numel(h_vals));
for j = 1:numel(Tau0_vals)
    tau0 = Tau0_vals(j);
    disp(['tau0 = ', num2str(tau0)]);
    guess = out(j,1);
    for i = 1:numel(h_vals)
        h = h_vals(i);
        params = poly_case(a0, alpha, tau0, R, h);
        params.marker = 's';
        params.pen = set_pen(...
            pens.lc{min(1, numel(pens.lc))}, ...
            pens.style{min(j, numel(pens.style))});
        %% exact
        mesh = set_full_mesh(N, params, 0);
        starter = @(sigma, params) starter_full(...
            @(problem, mesh)solver(problem, mesh, params), sigma, params, mesh);
        sigma = fit_sigma(starter, params, guess);
        %% save sigma to the output matrix
        out(j,i) = sigma;
        %% remeber last sigma as a guess
        guess = sigma;
    end
end