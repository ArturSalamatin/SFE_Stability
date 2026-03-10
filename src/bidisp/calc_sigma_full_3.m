function out = calc_sigma_full_3(...
    solver, N, alpha, Tau0_vals, a0, R, B_vals, pens)

zeta2 = zeros(size(Tau0_vals));
for k = 1:numel(zeta2)
    tau0_ = Tau0_vals(k);
    params = poly_case(a0, alpha, tau0_);
    zeta2(k) = z2(params);
end
h_vals = B_vals'*zeta2;

%% calculate sigma
out = -2 + alpha./(alpha+(1-alpha).*a0) ...
    + zeros(numel(B_vals), numel(Tau0_vals));
for i = 1:numel(B_vals)
        guess = out(i,1);
    for j = 1:numel(Tau0_vals)
        tau0 = Tau0_vals(j);
        disp(['tau0 = ', num2str(tau0)]);
        h = h_vals(i,j);
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
        out(i,j) = sigma;
        %% remeber last sigma as a guess
        guess = sigma;
        %% stop calculation for sufficiently large sigma
        if(guess > 0.7)
            out(i, (j+1):end) = guess;
            break;
        end
    end
end