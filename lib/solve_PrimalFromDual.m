function [R,szep] = solve_PrimalFromDual( Z )
% [R,szep] = solve_PrimalFromDual( Z )

% Decompose penalized matrix
[U,S,~] = svd(full(Z));
s = diag(S);

% Check SZEP condition,
szep = s(end)<1e-3 && (s(end)/s(end-1))<1e-3;
% Zero eigenvalue and margin to distinguish from numerical accuracy

if szep
  % Recovering the solution is trivial, just scale the eigenvector
  r = makenonhom( U(:,end) );
  R = reshape(r,3,3);
else
  error('No solution for non-tight case through nullspace yet')
end

end