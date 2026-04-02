function out = BC(params, eqN)
f = params.f;
%% BC at the left end
left = zeros(eqN,eqN);
left(1,1) = 1; % Phi(0) = 0
left(2,2) = 1; % Omega(0) = 0 /* = Psi(0)*/
% left(3,3) = 1; % Y(0) = 0
% left(4,4) = 1; % G(0)   = 1
%% BC at the right end
right = zeros(eqN,eqN);
% right(2,2) = 1; % Omega(1) = 0
right(3,1) = f; % f*Phi(1) + G(1) = 0
right(3,4) = 1; 
right(4,4) = 1; % G(1)   = 1
%% rhs for BC eqns
out.left = left;
out.right = right;
out.rhs = [0;0;0;1];
end