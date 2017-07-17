function [R,t,dstar,times,R0] = method_Olsson( correspondences )

tic
% Compute equivalent compressed quadratic form
q = compress_quadData( correspondences );
% normalize with the number of correspondences (mean error)
% q = (1/numel(correspondences))*q;
times.compress = toc;
% Solve the dual SDP problem
[dstar,Z,G,lam,time_cvx] = solve_dual_template( q, 'header_Olsson' );
times.cvx = time_cvx;
tic
% recover from primal SDP (dehomogeneize indexes for vec(R) and y)
[V,D] = eigs(G,2);
R0 = reshape(makenonhom(V([1:9,end],1)),3,3);
R = project_rotation( R0 );
t = makenonhom(V([10:12,end],1));
times.primalFromDual = toc;
times.solve = times.compress + times.cvx + times.primalFromDual;

end