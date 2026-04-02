function params = poly_case(a0, alpha, R, B)

%% set the packed bed
params.a0 = a0;
params.r = alpha; % dust volume fraction

params.a1 = 1.0;
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;

params.R = R;
params.B = B;
end
