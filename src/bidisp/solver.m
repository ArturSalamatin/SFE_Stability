function sol = solver(problem, mesh)
eqN = problem.eqN;
ids = problem.ids;
M = problem.M;
t = mesh.t;
I = mesh.I; % jump node is skipped
%% allocate memory
A = spalloc(M,M,M*2); % matrix
% A = full(A);
b = zeros(M, 1); % rhs
%% fill in the matrix
for i = I
    block_pos = (i-1)*eqN;
    step = (t(i+1) - t(i))/2.0; % (!)trapezoidal rule
    block = problem.block_matrix(t(i), t(i+1), i);
    Diag = problem.diag(i);
    % coef at y_i
    A(block_pos+ids, block_pos + ids) = ...
        block + Diag/step;
    % coef at y_(i+1)
    A(block_pos+ids, block_pos + ids + eqN) = ...
        block - Diag/step;
end
% BC -- boundary conditions
block_pos = (mesh.N-1)*eqN;
bc = problem.BC();
A(block_pos+ids, ids) = bc.left;
A(block_pos+ids, block_pos+ids) = bc.right;
b(block_pos+ids) = bc.rhs;
% JC -- jump condition
block_pos = (mesh.left.N-1)*eqN;
jc = problem.JC();
A(block_pos+ids, block_pos+ids) = jc.left;
A(block_pos+ids, block_pos+eqN+ids) = jc.right;
%% solution
sol.y = reshape(A\b, problem.eqN, mesh.N)';
sol.t = mesh.t;
% A = full(A);
end