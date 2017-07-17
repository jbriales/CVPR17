classdef Correspondence < handle & matlab.mixin.Heterogeneous
  %Correspondence Object with two 3D corresponding features
  %   See also Point2Point, Point2Line, Point2Plane.

%   Technical:
%   This is a *heterogeneous* object!
%   See also handle, matlab.mixin.Heterogeneous.

methods (Abstract)
  d2 = cost(this,T)
  % Compute cost for this correspondence with transformation T
  
  q = quad(this)
  % Return quadratic form for this correspondence (in terms of vec(T))
end

end