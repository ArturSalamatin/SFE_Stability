function out = Matrix(tau, z, params)
R = params.R;
h = params.z_step;
h2 = h*h;
alpha2 = (params.alpha)^2;

iter = params.iter;
x = params.x(:, iter);
dcdz = params.dcdz(:, iter);

m = numel(z);

I = 1:m;
J = I;
out = ...
    sparse(...
    [I,             I,        I(2:end)], ...
    [J,             J+m,      m+J(1:end-1)], ...
    [x.*dcdz*tau/2; g(x)+x/h; -x(2:end)/h], ...
    2*m, 2*m, 7*m) + ...
    sparse( ...
    m+[I,                I,            I(2:end),           I(1:end-1)], ...
    [J,                  J+m,          J(1:end-1),         J(2:end)], ...
    [2/h2+alpha2+R*dcdz; R*alpha2/tau*ones(m, 1); -1/h2*ones(m-1, 1); -(1/h2+R*dcdz(1:end-1))], ...
    2*m, 2*m, 7*m);

% inlet BC for psi
out(1, m+1) = 1; 
out(1, 1) = 0; 
% inlet BC for phi
out(m+1,1)   = 1;
out(m+1,2)   = 0;
out(m+1,m+1) = 0;
% inlet BC for phi
out(2*m,m)   = 1;
out(2*m,m-1) = 0;
out(2*m,2*m) = 0;

% A = full(out);

end