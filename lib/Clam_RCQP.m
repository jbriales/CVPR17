classdef Clam_RCQP
  %Clam_RCQP An object class for dual variables
  %   in the Lagrangian relaxation of the RCQP problem
  %   
  %   Detailed explanation goes here (TODO)
  
  properties
    unit_cols
    unit_rows
    ort_cols
    ort_rows
    det
    g
  end
  
  methods
    function this = Clam_RCQP(...
            lam_unit_cols,lam_unit_rows,...
            lam_ort_cols,lam_ort_rows,...
            lam_det,lam_g)

      if nargin==1
        % Vector-like input
        lamVec = lam_unit_cols;
        if ~isstruct(lamVec) && isvector(lamVec)
          assert(numel(lamVec)==22);
          lam_unit_cols = lamVec(1:3);
          lam_unit_rows = lamVec(4:6);
          lam_ort_cols  = lamVec(7:9);
          lam_ort_rows  = lamVec(10:12);
          lam_det       = reshape(lamVec(13:21),3,3);
          lam_g         = lamVec(22);
        else
          error('Wrong input for Clam_RCQP');
        end
      end
          
      % Assign values to fields
      this.unit_cols = lam_unit_cols;
      this.unit_rows = lam_unit_rows;
      this.ort_cols = lam_ort_cols;
      this.ort_rows = lam_ort_rows;
      this.det = lam_det;
      this.g = lam_g;
    end
    
    function maskedThis = mask(this, inputMask)
      if isa(inputMask,'Clam_RCQP')
        inputMask = inputMask.vec();
      end
      assert(islogical(inputMask));
      lamVec = this.vec();
      lamVec(inputMask) = 0;
      maskedThis = Clam_RCQP(lamVec);
    end
    
    function vOut = vec(this)
      vOut = [this.unit_cols;this.unit_rows;
              this.ort_cols;this.ort_rows;
              vec(this.det);
              this.g];
    end
    
    function disp(this)
      fields = fieldnames(this);
      for i=1:numel(fields)
        l = fields{i};
        fprintf([l,':\n']);
        if isvector(this.(l))
          disp(vec(this.(l))');
        else
          disp(this.(l));
        end
      end
    end
  end
  
  methods(Static)
    function lamMask = false()
      lamMask = Clam_RCQP( false(22,1) );
    end
    function lamMask = true()
      lamMask = Clam_RCQP( true(22,1) );
    end
    function lamMask = mask_cols()
      lam_unit_cols = true(3,1);
      lam_unit_rows = false(3,1);
      lam_ort_cols = true(3,1);
      lam_ort_rows = false(3,1);
      lam_det = false(3,3);
      lam_g = true;
      lamMask = Clam_RCQP(...
        lam_unit_cols,lam_unit_rows,...
        lam_ort_cols,lam_ort_rows,...
        lam_det,lam_g);
    end
    function lamMask = mask_rows()
      lam_unit_cols = false(3,1);
      lam_unit_rows = true(3,1);
      lam_ort_cols = false(3,1);
      lam_ort_rows = true(3,1);
      lam_det = false(3,3);
      lam_g = true;
      lamMask = Clam_RCQP(...
        lam_unit_cols,lam_unit_rows,...
        lam_ort_cols,lam_ort_rows,...
        lam_det,lam_g);
    end
    function lamMask = mask_colsAndRows()
      lam_unit_cols = true(3,1);
      lam_unit_rows = true(3,1);
      lam_ort_cols = true(3,1);
      lam_ort_rows = true(3,1);
      lam_det = false(3,3);
      lam_g = true;
      lamMask = Clam_RCQP(...
        lam_unit_cols,lam_unit_rows,...
        lam_ort_cols,lam_ort_rows,...
        lam_det,lam_g);
    end
    function lamMask = mask_all()
      lamMask = Clam_RCQP.true();
    end
  end
  
end