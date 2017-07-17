function [correspondences,T] = rand_registration( numObs, sd, R )
% [correspondences,T] = rand_registration( numObs, sd, R )
% Create a set of random correspondences in a sphere or radius R
% The *GT* transformation T used to generate the ideal problem is returned

%#ok<*NASGU>

if ~exist('numObs','var')
  numObs = 100;
end
if ~exist('sd','var')
  sd = 0.1;
%   sd = 0;
end
if ~exist('R','var')
  R = 10;
end
randMat = @() randn(3,numObs);

if isscalar(numObs)
  numObs = numObs * [1 1 1];
else assert(numel(numObs)==3)
end

% Generate random set of plane correspondences
p2p_correspondences = rand_pointToPoint(numObs(1),R);
p2l_correspondences = rand_pointToLine(numObs(2),R);
p2pl_correspondences = rand_pointToPlane(numObs(3),R);
% Concatenate all correspondences
correspondences = [p2p_correspondences, p2l_correspondences, p2pl_correspondences];
% Sanity check: cost with identity correspondence should be 0
% norm(correspondences.cost(Pose()))

% Generate random transformation
% This transformation convert measured points in frame {meas}
% into the coordinate frame of the model {model},
% so T = ^{model}T_{meas}
T = Pose.rand();

% Transform all sampled points into measured points (with inverse of GT)
allPoints = [correspondences.point];
transform(allPoints,inv(T))
% Sanity check: cost with GT transformation should be 0
% norm(correspondences.cost(T))

% Corrupt sampled points with noise
for i=1:numel(correspondences)
  % corrupt point coordinates with isotropic Gaussian
  corrupt(correspondences(i).point,sd);
end

end

function correspondences = rand_pointToPoint(m,R)

% preallocate objects
correspondences = repmat(Point2Point(),1,m);

for i=1:m % for each correspondence
  % Choose plane randomly such that it intersects the sphere
  % centered around zero with radius R.
  p1 = Point(rand_sphere(1,R));
  % The corresponding point is the same point
  p2 = copy(p1);
  
  % store new object
  correspondences(i) = Point2Point(p2,p1);
end
end

function correspondences = rand_pointToLine(m,R)

% preallocate objects
correspondences = repmat(Point2Line(),1,m);

for i=1:m % for each correspondence
  % Choose plane randomly such that it intersects the sphere
  % centered around zero with radius R.
  line = Line(rand_sphere(1,R),snormalize(randn(3,1)));
  
  % Sample a random point in the line
  % 1. Generate random 3D point inside the sphere
  rnd_point = Point(rand_sphere(1,R));
  % 2. Project point into line
  point = line.project(rnd_point);
  
  % store new object
  correspondences(i) = Point2Line(point,line);
end
end

function correspondences = rand_pointToPlane(m,R)

% preallocate objects
correspondences = repmat(Point2Plane(),1,m);

for i=1:m % for each correspondence
  % Choose plane randomly such that it intersects the sphere
  % centered around zero with radius R.
  plane = Plane(rand_sphere(1,R),snormalize(randn(3,1)));
  
  % Sample a random point in the plane
  % 1. Generate random 3D point inside the sphere
  rnd_point = Point(rand_sphere(1,R));
  % 2. Project point into plane
  point = plane.project(rnd_point);
  
  % store new object
  correspondences(i) = Point2Plane(point,plane);
end
end

function points = rand_sphere(m,R)

% directions (from normal distribution)
dir  = snormalize(randn(3,m));
% distances between 0 and R (uniform)
dist = R*rand(1,m);

% compose points
points = repmat(dist,3,1).*dir;

end