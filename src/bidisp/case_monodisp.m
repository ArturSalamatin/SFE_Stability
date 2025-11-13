function [params, sigma] = case_monodisp
params.a0 = 0.5;
params.a1 = 1.0;
params.r = 0.0; % dust volume fraction
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
params.R = 0.5;
params.h = 5;
params.t = 0.2;
params.a = sqrt(2*params.t);
params.z0 = z0(params);
params.z2 = z2(params);
params.dz2dt = dz2dt(params);
params.C1 = params.a./z2(params);
params.C2 = params.C1*params.a.*params.dz2dt;

sigma = -1.9072018;
end