function [R,t,times] = method_BnB( correspondences )

tic
[R,t] = BnB_solve(correspondences,1e-7);
times.solve = toc;

end