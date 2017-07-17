% Run this script to check/set dependencies in the path

% Add geometry toolkit
addpath(fullfile(pwd,'geometry'))
% Add methods for convex relaxation of rotation
addpath(fullfile(pwd,'lib'))
% Add solvers (specific to registration)
addpath(genpath(fullfile(pwd,'solvers')))
% Add experimental real data
addpath(fullfile(pwd,'data'))

% Add toolbox of algebra functions
addpath(fullfile(pwd,'mMath'))

% CVX: modelling tool for convex optimization
if ~exist('cvx_begin')
  error('CVX not available, install toolbox');
end

% SeDuMi: solver for SDP (called by BnB solver)
if ~exist('sedumi')
  warning('sedumi not available as an independente toolbox, required by BnB solver');
end

disp('SETUP: All dependencies added correctly')