classdef Line < Feature
  %Line Object with 3D line data fields (point and direction)
  %   See also Feature.
  
  %#ok<*PROPLC>
  
  properties
    x % 3D point in the line (center?)
    v % 3D direction (in S2)
  end
  
  methods
    function this = Line( point,direction )
      % line = Line( point,direction )
      %   Constructor from components
      
      if nargin == 0
        % default values: vertical line at Z=0
        point  = [0 0 0]';
        direction = [0 0 1]';
      end

      this.x = point;
      this.v = direction;
    end

    function point_proj = project(this,point)
      % x_proj = project( LINE, POINT )
      % Project POINT into LINE by finding its closest point
      
      assert( isa(point,'Point'), 'Invalid input point' );
      
      % A point x is projected to a line (p,v) by p + vv'*(x-p)
      parDist = dot(this.v,point.x-this.x,1);
      x_proj = this.x + this.v.*(ones(3,1)*parDist);
      % Sanity check: sample points must belong to the lines: v'*(x-p)=0
      % norm(cross(this.v,x_proj-this.x,1))
      
      point_proj = Point(x_proj);
    end
    
    function d2 = distSq(this,point)
      % d^2 = distSq(LINE,POINT)
      % Compute squared distance from POINT to LINE
      
      closestDist = norm( (eye(3)-this.v*this.v')*(point.x-this.x) );
      d2 = closestDist^2;
    end
    
    function h = plot(this,varargin)
      color = 'k';
      hP = plot3(this.x(1),this.x(2),this.x(3),[color,'x']);
      % Set line (1m for now)
      L = 1;
      lineData = [ this.x, this.x+L*this.v ]';
      c_foo = num2cell(lineData,1);
      hL = line(c_foo{:},'Color',color);
      
      h = [hP;hL]; % concatenate handles
    end
  end
end