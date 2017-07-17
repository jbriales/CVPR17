% Main script
% Example to solve different multi-modal registration problems
% with the method proposed in [1], and the previous approaches in [2,3]
% 
% Related work:
% [1] Briales, J., & Gonzalez-Jimenez, J. "Convex Global 3D Registration with Lagrangian Duality." In CVPR 2017
% [2] Olsson, C., & Eriksson, A. "Solving quadratically constrained geometrical problems using lagrangian duality." In ICPR 2008.
% [3] Olsson, C., Kahl, F., & Oskarsson, M. "Branch-and-Bound Methods for Euclidean Registration Problems." In TPAMI 2009.

clear all

% choose the problem to solve
problemType = 'random';
% problemType = 'SpaceStation';
% problemType = 'RubikCube';

% Problem generation
% -------------------------------------------------------------------------
switch problemType
  case 'random'
    % create random problem
    % number of correspondences: [nr points, nr lines, nr planes]
    problem.m = [2 4 7]; 
    % noise in the correspondences
    problem.noise = 0.1;
    % size of the random scene
    scene_radius = 10;
    [correspondences,gt_T] = rand_registration( problem.m, problem.noise, scene_radius );
    
  case 'SpaceStation'
    % registration data from Space Station [2,3]
    [c_p,c_l,c_pl] = correspondences_SpaceStation( );
    correspondences = [c_p,c_l,c_pl];
    
  case 'RubikCube'
    % registration data from Rubik cube [2,3]
    [c_p,c_l,c_pl] = correspondences_RubikCube( );
    correspondences = [c_p,c_l,c_pl];
    
  otherwise
    error('Unknown problem type')
end

% Model the problem as a compressed quadratic form in R only
% -------------------------------------------------------------------------
% Compute equivalent compressed quadratic form
q = compress_quadData( correspondences );
sv = svd(q.Q_);
% Sanity check: Sum of costs and compressed quad form must be equivalent
% abs( q.eval(vec(gt_T)) - sum(cost(correspondences,gt_T)) )

% Marginalize the quadratic function wrt translation using Schur complement
t_idxs = 10:12;
[q_margin, A] = marginalize(q,t_idxs);

% Solve the problem with our method
% ------------------------------------------------------------------------
[R,t,dstar,times] = method_RCQP( correspondences, 'header_all' );
T = Pose(t,R);

f = q.eval(vec(T));
gap = (f-dstar)/dstar;
fprintf('Optimality gap is f^star-d^star=%E\n',gap);

% Solution from Olsson
% -------------------------------------------------------------------------
[R,t,d,times] = method_Olsson( correspondences );
T = Pose(t,R);

f = q.eval(vec(T));
gap2 = (f-d)/d;
fprintf('Optimality gap is f^star-d^star=%E\n',gap2);

% Solution with BnB
% -------------------------------------------------------------------------
tic
[R_BnB,t_BnB] = BnB_solve(correspondences,1e-10);
time_BnB = toc;
T_BnB = Pose(t_BnB,R_BnB);
f_BnB = q.eval(vec(T_BnB));
gap3 = (f_BnB-dstar)/dstar;
fprintf('Optimality gap is f^star-d^star=%E\n',gap3);