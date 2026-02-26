function params = poly_case(a0, alpha, t, R, h)
if(nargin < 5)
    h = 0;
end
if(nargin < 4)
    R = 0;
end
%% set the packed bed
params.a0 = a0;
params.r = alpha; % dust volume fraction

params.a1 = 1.0;
params.g1 = (1-params.r)/params.a1;
params.g0 = params.g1 + params.r/params.a0;
%% set time
if(nargin == 2)
    t = (a0*a0/2+0.5)/2;
end
if(t >= 0.5)
    error("Time moment is too big!");
end
if(t < a0*a0/2)
    error("Time moment is too small!");
end
params.t = t;
params.R = R;
params.h = h;
%% set dependent vars
params.a = sqrt(2*t);
params.z0 = z0(params);
params.z2 = z2(params);
params.xi0 = params.z0/params.z2;
params.dz2dt = dz2dt(params);
params.dz0dt = dz0dt(params);
params.C1 = params.a./params.z2;
params.C2 = params.C1*params.a.*params.dz2dt;
end
