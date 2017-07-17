classdef Quadratic < handle & matlab.mixin.Copyable
%Quadratic A class to treat quadratic functions as an object
% Access usual fields for a quadratic function as a struct:
%   Q - quadratic part
%   b - linear part
%   c - constant term
% Also available:
%   Q_- homogeneized version of the quadratic form:
%       | Q  b |
%       | b' c |
% 
% Constructor:
%   q = Quadratic(Q,b,c)
%   q = Quadratic(Q_)
% 
% Convenient methods:
%   eval(q,x)
  
  properties
    Q
    b
    c
  end
  properties (Dependent)
    Q_
  end
  
  methods
    function this = Quadratic(varargin)
      % Quadratic(Q,b,c)
      % Quadratic(Q_)
      % Terms for a quadratic function with matrix
      %  [ Q  b
      %    b' c ]
      % NOTE: b is the halved value of the linear term!
      
      if nargin == 1
        % Quadratic(Q_)
        this.Q_ = varargin{1};
      elseif nargin == 3
        % Quadratic(Q,b,c)
        this.Q = varargin{1};
        this.b = varargin{2};
        this.c = varargin{3};
      else
        error('Wrong number of arguments');
      end
      % Fill in with zeros by default
      if isempty(this.Q); this.Q = zeros(9); end
      if isempty(this.b); this.b = zeros(9,1); end
      if isempty(this.c); this.c = 0; end
      
      assert(issymmetric(this.Q_),'Assert input is symmetric')
    end
    % Interface the complete matrix
    function Q_ = get.Q_(this)
      Q_ = [this.Q  this.b
            this.b' this.c];
    end
    function set.Q_(this,Q_)
      assert(issymmetric(Q_),'The quadratic matrix is not symmetric');
      this.Q = Q_(1:end-1,1:end-1);
      this.b = Q_(1:end-1,end);
      this.c = Q_(end,end);
    end
    % Convenient test function
    function hom = ishom(this)
      if all(this.b==0) && this.c==0
        hom = true;
      else
        hom = false;
      end
    end
    
    function f = eval(this, x)
      % Evaluate the current quadratic function at the vector x
      % q(x) = x'*Q*x + 2*b'*x + c = hom(x)'*Q_*hom(x)
      assert(isnumeric(x) && isvector(x),'Input should be a vector')
      f = makehom(x)'*this.Q_*makehom(x);
    end
    
    function [q_margin, A] = marginalize(q,i)
      % [q_margin, A] = marginalize(q,MARGINALIZED_IDXS)
      % Compute the marginalization of a quadratic form
      % wrt the variables y in positions MARGINALIZED_IDXS.
      % This consists of finding the critical point wrt those variables,
      % find the optimal value in terms of the rest of unknowns as
      %   y_optim = A*x
      % and obtain the quadratic form q(x) wrt the remaining variables x.
      %
      % The results here are obtained from directly deriving the quadratic
      % function (simpler if homogeneized) wrt the desired variables.
      % The results are equivalent to applying the Schur complement.
      
      % Get matrix data
      M = q.Q_;
      dim = size(M,1);
      not_i = setdiff(1:dim,i);
      
      % Compute marginalized cost function (a new quadratic form)
      M_margin = M(not_i,not_i) - M(not_i,i)*pinv(M(i,i))*M(i,not_i);
      q_margin = Quadratic(symmetrize(M_margin));
      
      if nargout == 2
        % Compute linear map from non-marginalized x to marginalized y
        s = svd(M(i,i));
        if s(end) < 1e-2
          warning('Marginalization possibly degenerate: sv(end)=%E',s(end))
        end
        A = -pinv(M(i,i))*M(i,not_i);
      end
    end
    
    function q = plus(a,b)
      % q = plus(a,b)
      % q = a+b
      % Compute the sum of the quadratic matrices underlying a and b
      new_Q = symmetrize(a.Q_ + b.Q_); % symmetrize for num stability
      q = Quadratic( new_Q );
    end
    function r = mtimes(a,q)
      assert( isnumeric(a) && isscalar(a), 'Invalid scalar-quad product' )
      r = Quadratic( a*q.Q_ );
    end
        
    function disp(this)
      if nargin > 1
        for i=1:numel(this)
          disp(this(i));
        end
        return
      end
      disp(this.Q_)
    end
  end
  
  
end