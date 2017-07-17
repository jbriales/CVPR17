function [d,Z,G,lam,time] = solve_dual_template( q, header )
% [d,Z,G,lam,time] = solve_dual_template( q )
% [d,Z,G,lam,time] = solve_dual_template( q, header )
% 
% Generic solver for the dual of a Rotation-Constrained Quadratic Program
% of the form
%   min_{R\inSO(3)} r'Qr+b'r+c (RCQP)
% using only quadratic constraints (leads to QCQP primal formulation)
% 
% The different relaxation options can be chosen by using
% different header scripts when defining the dual variables:
% That is, for each approach, dual variables are defined or set to zero

%#ok<*VUNUS>

if nargin < 2
  header = 'header_all';
end

% Set CVX options
cvx_precision('best') % TODO: Choose as option from outside?
cvx_solver sedumi % TODO: Choose as option from outside?
% cvx_solver SDPT3
cvx_quiet(false) % necessary false to dump results
dumpfile = sprintf('temp_cvxDump_%s',num2str(randi(1e10)));
cvx_solver_settings( 'dumpfile', dumpfile ) % dump statistics for parsing

tic
cvx_begin SDP % Use CVX's SDP mode
%% Define all dual variables and penalize
% Use a script to have flexibility in the definition here
% and keep compatibility with the CVX environment
run(header)
% Slack variable (for Schur complement or homogeneizing quad form)
variable lam_g
% Stack variables into structure
lam = Clam_RCQP( lam_unit_cols, lam_unit_rows,...
                 lam_ort_cols, lam_ort_rows,...
                 lam_det, lam_g );

% Compute penalization matrix
if ~exist('P','var') % this allows defining custom P inside header
  P = build_penalizationMat( lam );
else
  P = P + build_penalizationMat( lam );
end
% Build slack matrix: -P + Z = Q → Z = Q + P
% Penalize homogeneous quadratic matrix (PSD slack matrix)
if all(size(q.Q_)==[10 10])
  % Marginalized case: [vec(R) y]
  Z = q.Q_ + P;
elseif all(size(q.Q_)==[13 13])
  % Non-marginalized case: [vec(R) t y]
  % Pad with zero blocks in penalization
  P_padded = [ P(1:9,1:9) zeros(9,3) P(1:9,end)
               zeros(3,9) zeros(3,3) zeros(3,1)
               P(end,1:9) zeros(1,3) P(end,end) ];
	Z = q.Q_ + P_padded;
else
  error('Wrong dimensions of the quadratic function')
end
% Set SDP objective and constraints
Z = symmetrize(Z);
dual variable G
G : Z >= 0; % Slack matrix must be positive semidefinite
d = lam_g;
maximize(d)
cvx_end
fprintf('CVX-Matlab time: %E\n',toc);

% Reload CVX results to output object
lam = Clam_RCQP( lam_unit_cols, lam_unit_rows,...
                 lam_ort_cols, lam_ort_rows,...
                 lam_det, lam_g );
               
% Read execution time
CPU_total = readDump_time(dumpfile);
time = CPU_total;
delete([dumpfile,'.mat'])
% For clearness, unset dumping as solver option
cvx_solver_settings -clear dumpfile

end

% Auxiliar function to parse time from SDPT3 log
function CPU_total = readDump_time( dumpFile )
% readDump_time( outputStr )
% Read the CPU time from screen output
% Currently assumes SDPT3 format

load(dumpFile, 'output')
% Check kind of solver
solverName = cvx_solver;
switch solverName
  case 'SeDuMi'
    % Find detailed timing and split char
    m = regexp(output, 'Detailed timing', 'split');
    m = regexp(m{2}, 'Post', 'split');
    m = regexp(m{2}, 'Max-norms', 'split');
    str = m{1}(2:end-1); % Skip newlines
    times = sscanf(str,'%E');
    t.Pre = times(1);
    t.IPM = times(2);
    t.Post= times(3);
    CPU_total = sum(times);
  case 'SDPT3'
    % Find nearby lines
    m = regexp(output, 'Total CPU time .* \n', 'match');
    str = m{1};
    % Regexp for numbers in the two lines
    times = regexp(str,'\d.\d*','match');
    % Convert char numbers for numeric values
    times = cellfun(@str2num,times);
    
    CPU_total = times(1);
%     CPU_it = times(2);
  otherwise
    error('Unrecognized solver');
end

end
