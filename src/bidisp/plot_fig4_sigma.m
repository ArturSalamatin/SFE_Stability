clc
clear all
% close all

global sigma_fig
sigma_fig = 6;
%% packed bed params
a0 = 0.2;
alpha = [0, 0.1, 0.3, 0.5, 0.7];
marker = ['o','s','d','*','x'];
%% time
t = [0.2, 0.03];
%% perturbation parameters
h = [0.1, 0.5, 1, 2, 3, 5, 10];% linspace(1,505,11);
R = linspace(0.01, 3, 85);
col = {'m'};
style = {'-'};
%% set mesh
N = 300;
q = 1.06;
for t_id = 1:numel(t)
    sigma_fig = sigma_fig+1;
    my_figure(sigma_fig)
    hold on
    axis([0 3 -3.5 0])
    % set curves to calculate and plot
    if(t_id == 1)
        alpha_id = 1:5;
    elseif(t_id == 2)
        alpha_id = [1,5];
    end
    for J = alpha_id
        %% set the packed bed
        params = poly_case(a0, alpha(J), t(t_id));
        params.marker = marker(J);
        for i = 1:numel(N)
            %         mesh = quasiuniform_mesh(N, params);
            mesh = log_mesh(q, N(i), params);
            %% choose starter
            starter = @(sigma, params) starter_Y(...
                sigma, params, mesh);
            pen = set_pen(col, style);
            out = calc_sigma(h, R, starter, pen, params);
        end
    end
end


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
