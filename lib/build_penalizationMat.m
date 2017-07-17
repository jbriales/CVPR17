function Qpen_ = build_penalizationMat( lam )
% Qpen_ = build_penalizationMat( lam )
% 
% Creates the *homogeneous* penalization matrix for LMI in SDP:
%   Z = Q0_ + Qpen_(lam) >= 0
% 
% The penalization matrix depends *only* on the dual variables lam
% It is independent of the problem data or the primal variables
% 
% The expressions are mat-like (more compact), see the report for details 

% Assert the values of lam are given struct-like
if numel(lam)==22 || ~isstruct(lam) && ~isobject(lam)
  % Convert from vec-format to struct-format
  lam = Clam_RCQP(lam);
end

% Quadratic term
Lam_cols = diag(lam.unit_cols) + 0.5*offsym(lam.ort_cols);
Lam_rows = diag(lam.unit_rows) + 0.5*offsym(lam.ort_rows);
skew_det = cell(1,3);
for k=1:3
  skew_det{k} = skew(lam.det(:,k));
end
Q_pen = -kron(Lam_cols,eye(3)) - kron(eye(3),Lam_rows) + 0.5*skew( skew_det );

% Linear term
b_pen = - 0.5*vec(lam.det);

% Independent term
c_pen = sum(lam.unit_cols) + sum(lam.unit_rows) - lam.g;

% Complete homogeneous matrix
Qpen_ = [Q_pen  b_pen;
         b_pen' c_pen];

end