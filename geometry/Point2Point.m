classdef Point2Point < Point2X
  %Point2Point Correspondence between 3D point and 3D point
  %   See also Correspondence.
  
  %#ok<*PROPLC>
  
  properties
    % rank? 3
    % + inherited from Point2X
    format = 'r*';
  end
  
  methods
    function this = Point2Point( point,modelPoint )
      % correspondence = Point2Point( POINT,MODEL_POINT )
      %   Constructor from components
      
      if nargin == 0
        % default values
        point = Point();
        modelPoint = Point();
      end

      this.point = point;
      this.model = modelPoint;
    end
    
    function q = quad(this)
      % Return quadratic form for this correspondence (in terms of vec(T))
      
      % Define point concatenation matrix
      Xy = [ kron(makehom(this.point.x)',eye(3)), -this.model.x ];
      % Define central matrix for point-to-point
      C  = eye(3);
      
      % Create homogeneized quadratic in [vec(R),t,1]=[vec(T),1]
      Q_ = symmetrize(Xy'*C*Xy); % symmetrize for num stability
      q = Quadratic(Q_);
    end
    
  end
end
