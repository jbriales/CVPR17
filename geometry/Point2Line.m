classdef Point2Line < Point2X
  %Point2Line Correspondence between 3D point and 3D line
  %   See also Correspondence.
  
  %#ok<*PROPLC>
  
  properties
    % inherited from Point2X
    format = 'go';
  end
  
  methods
    function this = Point2Line( point,modelLine )
      % correspondence = Point2Line( POINT,MODEL_LINE )
      %   Constructor from components
      
      if nargin == 0
        % default values
        point = Point();
        modelLine  = Line();
      end

      this.point = point;
      this.model  = modelLine;
    end
    
    function q = quad(this)
      % Return quadratic form for this correspondence (in terms of vec(T))
      
      % Define point concatenation matrix
      Xy = [ kron(makehom(this.point.x)',eye(3)), -this.model.x ];
      % Define central matrix for point-to-model
      C  = eye(3) - this.model.v*this.model.v';
      
      % Create homogeneized quadratic in [vec(R),t,1]=[vec(T),1]
      Q_ = symmetrize(Xy'*C*Xy); % symmetrize for num stability
      q = Quadratic(Q_);
    end
    
  end
end
