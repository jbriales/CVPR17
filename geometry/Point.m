classdef Point < Feature
  %Point Object with 3D point data fields (point only)
  %   See also Feature.
  
  %#ok<*PROPLC>
  
  properties
    x % 3D point
  end
  
  methods
    function this = Point( point )
      % plane = Point( point )
      %   Constructor from components
      
      if nargin == 0
        % default values: point at origin
        point  = [0 0 0]';
      end

      this.x = point;
    end
    
    function d2 = distSq(this,point)
      % d^2 = distSq(PLANE,POINT)
      % Compute squared distance from POINT to PLANE
      
      d2 = norm(this.x-point.x)^2;
    end
    
    % TODO: This should be available in all feature types? Then abstract.
    function new_point = corrupt(this,sd)
    % corrupt(point,noise)
    % new_feature = corrupt(point,noise)
    % Corrupt this feature with noise
    % TODO: For now Gaussian, but should this be generic function?
    
    % generate Gaussian random noise 3D vector
    new_point = Point(this.x+sd*randn(3,1));
    if nargout==0
      copyfields(this,new_point);
    end
    end
  
    function h = plot(this,varargin)
      h = plot3(this.x(1),this.x(2),this.x(3),varargin{:});
    end
    
  end
end
