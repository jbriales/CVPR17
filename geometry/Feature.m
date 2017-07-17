classdef Feature < handle & matlab.mixin.Copyable
  %Plane Object with 3D feature
  %   See also Point, Line, Plane.

%   Technical:
%   This is a *handle* *copyable* object, that means that the object
%   created by the constructor is a handle (pointer) to the data.
%   If you want to get an *independent* copy of the object (clone)
%   use the method *copy* in the object.
%   Intuitively, the pose of an object is inherent to this
%   and should be assigned to a different object, but another object
%   could have the same value (copy).
%   See also handle, matlab.mixin.Copyable.

methods
  function new_feature = transform(this,T)
    % transform(this,T)
    % new_feature = transform(this,T)
    % Transform this feature with given T
    % If nargout==0 it performs the transformation in-place (handle class)
    
    if numel(this)>1
      % handle array input
      if nargout==0
        for i=1:numel(this)
          transform(this(i),T);
        end
      else
        new_feature = this;
        for i=1:numel(this)
          new_feature(i) = transform(this(i),T);
        end
      end
      return
    end

    % use overloaded mtimes in transformation as composition operation
    new_feature = T*this; 
    if nargout==0
      % copy output data into the fields of the handle input
      copyfields(this,new_feature);
    end
  end
end
  
end