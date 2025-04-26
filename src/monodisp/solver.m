function sol = solver(problem, mesh)
eqN = problem.eqN;
ids = problem.ids;
M = problem.M;
t = mesh.t;
I = mesh.I;
%% allocate memory
Diag = sparse(1:eqN, 1:eqN, ones(1,eqN), eqN, eqN, eqN);
A = spalloc(M,M,M*2);
%% fill in the matrix
block_pos = 0;
for i = I
    step = (t(i+1) - t(i))/2; % (!)trapezoidal rule
    block = problem.block_matrix(t(i), t(i+1));
    % coef at y_i
    A(block_pos+ids, block_pos + ids) = ...
        block + Diag/step;
    % coef at y_(i+1)
    A(block_pos+ids, block_pos + ids + eqN) = ...
        block - Diag/step;
    % move to the next iteration
    block_pos = block_pos+eqN;
end
% BC
bc = problem.BC();
A(block_pos+ids, ids) = bc.left;
A(block_pos+ids, block_pos+ids) = bc.right;
%% rhs
b = zeros(M, 1);
b(block_pos+ids) = bc.rhs;
%% solution
sol.y = reshape(A\b, problem.eqN, mesh.N)';
sol.t = mesh.t;
end