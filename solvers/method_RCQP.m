function [R,t,dstar,times,R0] = method_RCQP( correspondences, header )

tic
% Compute equivalent compressed quadratic form
q = compress_quadData( correspondences );
% normalize with the number of correspondences (mean error)
% q = (1/numel(correspondences))*q;
times.compress = toc;
tic
% Marginalize the quadratic function wrt translation using Schur complement
[q_margin, A] = marginalize(q,10:12);
times.margin = toc;
% Solve the dual SDP problem
[dstar,Z,G,lam,time_cvx] = solve_dual_template( q_margin, header );
times.cvx = time_cvx;
tic
% Solve primal solution from dual via complementary slackness
R0 = solve_PrimalFromDual(Z);
% get chordal projection (closest rotation)
R = project_rotation( R0 );
% Recover optimal marginalized translation component
t = A*makehom(vec(R0));
% finish timer
times.primalFromDual = toc;
times.solve = times.compress + times.margin + times.cvx + times.primalFromDual;

end