%% TO BE CALLED FROM solve_dual_template function

% Normality constraints (unit norm)
variables lam_unit_cols(3,1) lam_unit_rows(3,1)
% Orthogonality constraints
variables lam_ort_cols(3,1) lam_ort_rows(3,1)
% Determinant constraints
variable lam_det(3,3) % Stack lam_det_ijk as column vectors in single mat