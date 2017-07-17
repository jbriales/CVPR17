classdef Point2X < Correspondence
  %Point2X Abstract correspondence between 3D point and 3D primitive
  %   See also Correspondence, Point2Point, Point2Line, Point2Plane.
  
  %#ok<*PROPLC>
  
  properties
    point % always a Point feature
    model % a Feature of the types Point, Line or Plane
  end
  
  methods (Sealed)
    function d2 = cost(this,T)
      % d2 = cost(CORRESPONDENCE,TRANSFORMATION)
      % Compute cost function for this correspondence:
      %   dist(T(point),point)
      
      if numel(this) > 1
        % handle array input
        d2 = NaN(size(this));
        for i=1:numel(this)
          d2(i) = cost(this(i),T);
        end
        return
      end
      
      d2 = distSq(this.model,T*this.point);
    end    
  end
end
