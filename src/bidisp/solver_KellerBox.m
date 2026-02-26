function sol = solver_KellerBox(problem, mesh, params)
%% ODE solver using Keller box scheme
eqN = problem.eqN; % nmbr of equations
ids = problem.ids; % ids within a block, 1:eqN
M = problem.M; % nmbr of discrete unknows
xi = mesh.xi;
%% allocate memory
A = spalloc(M,M,M*eqN); % matrix
b = zeros(M, 1); % rhs
% fill in the matrix
%% fill regular segments
segmIds = mesh.segmIds; % jump node is skipped
for i = segmIds
    block_pos = (i-1)*eqN;
    step = (xi(i+1) - xi(i))/2.0; % (!)centered Euler scheme, deHoog1976
    block = problem.block_matrix(i);
    Diag = problem.diag(i);
    % coef at y_i
    A(block_pos+ids, block_pos + ids) = ...
        block + Diag/step;
    % coef at y_(i+1)
    A(block_pos+ids, block_pos + ids + eqN) = ...
        block - Diag/step;
end
%% jump conditions, if any
jumpIds = mesh.jumpIds; % jump segments, may be empty
for i = jumpIds
    block_pos = (i-1)*eqN;
    jc = problem.JC();
    A(block_pos+ids, block_pos+ids) = jc.left;
    A(block_pos+ids, block_pos+eqN+ids) = jc.right;
    b(block_pos+ids) = jc.rhs;
end
%% boundary conditions
% BC -- boundary conditions
N = mesh.N;
block_pos = (N-1)*eqN;
bc = problem.BC_Keller();
A(block_pos+ids, ids) = bc.left;
A(block_pos+ids, block_pos+ids) = bc.right;
b(block_pos+ids) = bc.rhs;
%% solution
sol.id = mesh.left.N;
sol.y = reshape(A\b, problem.eqN, mesh.N)';
sol.t = mesh.xi;
sol.x = mesh.x;
sol.mesh = mesh;
%% scale as if X(jump_node) = 1
factor = sol.y(sol.id+1,2);
sol.y = sol.y/factor;
end