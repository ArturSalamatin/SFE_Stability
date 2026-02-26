function params = set_params(a0, alpha, R, B, H)
% size of dust particles
params.a0 = a0;
% volume fraction of dust particles
params.r = alpha; % dust volume fraction
% size of large particles
params.a1 = 1.0;
% SSA of large particles
params.g1 = (1-params.r)/params.a1;
% overall SSA
params.g0 = params.g1 + params.r/params.a0;
% ratio of two spatial scales
params.B = B;
% dimensionless vessel height
params.H = H;
% exponent for viscosity
params.R = R;
end

