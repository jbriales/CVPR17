%% TO BE CALLED FROM solve_dual_template function

% Orthonormal columns constraints
variables lam_unit_cols(3,1) lam_ort_cols(3,1)
% Set the rest of multipliers to zero (non-existent)
lam_unit_rows = zeros(3,1);
lam_ort_rows  = zeros(3,1);
lam_det       = zeros(3,3);