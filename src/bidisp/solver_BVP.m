function sol = solver_BVP(problem, xmesh, params)

opts = bvpset(...
    'RelTol', 1e-7...
    ,'AbsTol', 1e-7); 
%'SingularTerm', problem.BVP_S());

guess = solinit(problem, xmesh, params);

out = bvp5c(...
    @(x,y)problem.BVP_f(x,y), ...
    @(ya,yb)problem.BVP_bc(ya,yb), ...
    guess, opts);
% u == out.x == a - x
% figure(700)
% plot(params.a - out.x, out.y(2,:)/out.y(2,end), 'b-')
% plot(params.a - out.x, out.y(1,:)/out.y(2,end), 'b--')
sol.x = params.a - out.x;
sol.y = out.y'/out.y(2,end);
sol.t = z_of_x(sol.x, params)/params.z2;
sol.id = numel(sol.x);

end

function vals = pff(u,xmesh, vals, params)
x = min(xmesh(1),max(params.a - u, xmesh(end)));
vals = interp1(xmesh, vals, x);

if(any(isnan(vals)))
    warning('NaN interp.');
end

end

function sol_out = solinit(problem, xmesh, params)
a = params.a;
sol = solver_RK(problem, xmesh, params);

xmesh = [sol.x];
vals = [sol.y];

umesh = a - xmesh;

sol_out = bvpinit(umesh, @(u)pff(u, xmesh, vals, params));

% figure(700)
% hold on
% % plot(sol.x, sol.y, 'r')
% plot(a - umesh, sol_out.y(2,:)/sol_out.y(2,end), 'k-')
% plot(a - umesh, sol_out.y(1,:)/sol_out.y(2,end), 'k--')

end
