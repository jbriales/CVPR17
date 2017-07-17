% function BnB_setData( correspondences )
% function [B,k,m,mpi,ml,mp] = BnB_setData( correspondences )
% Better use as script

global m;
global B;
global k;
global mpi;
global ml;
global mp;

%% Convert correspondences into Olsson's format
% m - scalar, rest - 3xm matrices
% Point-to-plane
idxs = arrayfun(@(x)isa(x,'Point2Plane'),correspondences);
mpi = sum(idxs);
if mpi > 0
  allPoints = [correspondences(idxs).point];
  xpi = [allPoints.x];
  allPlanes = [correspondences(idxs).model];
  ypi = [allPlanes.x];
  n   = [allPlanes.n];
end
% Point-to-line
idxs = arrayfun(@(x)isa(x,'Point2Line'),correspondences);
ml = sum(idxs);
if ml > 0
  allPoints = [correspondences(idxs).point];
  xl = [allPoints.x];
  allLines  = [correspondences(idxs).model];
  yl = [allLines.x];
  v  = [allLines.v];
end
% Point-to-point
idxs = arrayfun(@(x)isa(x,'Point2Point'),correspondences);
mp = sum(idxs);
if mp > 0
  allPoints = [correspondences(idxs).point];
  xp = [allPoints.x];
  allModelPoints = [correspondences(idxs).model];
  yp = [allModelPoints.x];
end

%% Apply code from Olsson
% Read variables:
% Point-to-plane:
%   mpi, xpi, ypi, n
% Point-to-line:
%   ml, xl, yl, v
% Point-to-point:
%   mp, xp, yp

% Ber�kna matriser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point-to-point central matrix
N = zeros(3,3);
if mp > 0
   N = N + mp*eye(3);
end

% Point-to-line central matrix
if ml > 0
   V = cell(1,ml);
   for i = 1:ml;
      V{i} = eye(3) - v(:,i)*v(:,i)';
      N = N + V{i}'*V{i};
   end
end

% Point-to-plane central matrix
if mpi > 0
   for i = 1:mpi
      N = N + n(:,i)*n(:,i)';
   end
end

if rank(N) < 3;
    disp('N ej inverterbar');
end

% Ber�kna m�lfunktion (punkt-plan delen) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mpi > 0
   A = cell(1,mpi);
   for i = 1:mpi
      a = n(:,i);
      b = xpi(:,i);
      A{i} = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
	      a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
	      
	      a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
	      a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
	      
	      a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
	      -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		
	      a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
	      a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
   end
   
   
   A1 = cell(1,mpi);
   for i = 1:mpi
      A1{i} = zeros(4,4);
      slask = zeros(4,4);
      for j = 1:mp
	 a = (n(:,i)'*N^-1)' ;
	 b = xp(:,j);
	 slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		  
		  a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		  a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		  
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		  -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		  
		  a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		  a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	 A1{i} = A1{i} + slask;
      end
   end
   
   k1 = cell(1,mpi);
   for i = 1:mpi
      slask = 0;
      for j = 1:mp
	 slask = slask + n(:,i)'*N^-1*yp(:,j);
      end
      k1{i} = slask;
   end
   
   A2 = cell(1,mpi);
   for i = 1:mpi 
      A2{i} = zeros(4,4);
      slask = zeros(4,4);
      for j = 1:ml
	 a = (n(:,i)'*N^-1*V{j}'*V{j})' ;
	 b = xl(:,j);
	 slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		  
		  a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		  a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		  
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		  -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		  
		  a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		  a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	 A2{i} = A2{i} + slask;
      end
   end
   
   k2 = cell(1,mpi);
   for i = 1:mpi
      slask = 0;
      for j = 1:ml
	 slask = slask + n(:,i)'*N^-1*V{j}'*V{j}*yl(:,j);
      end
      k2{i} = slask;
   end
   
   A3 = cell(1,mpi);
   for i = 1:mpi
      slask = zeros(4,4); 
      A3{i} = zeros(4,4);
      for j = 1:mpi
	 a = (n(:,i)'*N^-1*n(:,j)*n(:,j'))' ;
	 b = xpi(:,j);
	 slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		  
		  a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		  a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		  
		  a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		  -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		  
		  a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		  a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	 A3{i} = A3{i} + slask;
      end
   end

   k3 = cell(1,mpi);
   for i = 1:mpi
      slask = 0;
      for j = 1:mpi
	 slask = slask + n(:,i)'*N^-1*n(:,j)*n(:,j)'*ypi(:,j);
      end
      k3{i} = slask;
   end
   
   B = cell(1,mpi);
   k = cell(1,mpi);
   for i = 1:mpi
      B{i} = A{i}-A1{i}-A2{i}-A3{i};
      k{i} = k1{i} + k2{i} + k3{i} - n(:,i)'*ypi(:,i);
   end
end

   
% Ber�kna m�lfunktion (punkt-linje delen) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ml > 0;
   A = cell(1,ml);
   for ii = 1:3
      for i = 1:ml
	 vik = V{i}(ii,:);
	 a = vik';
	 b = xl(:,i);
	 A{i} = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		 a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		 
		 a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		 a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		 
		 a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		 -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		 
		 a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		 a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
      end
      
      
      A1 = cell(1,ml);
      for i = 1:ml
	 vik = V{i}(ii,:);
	 A1{i} = zeros(4,4);
	 slask = zeros(4,4);
	 for j = 1:mp
	    a = (vik*N^-1)' ;
	    b = xp(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A1{i} = A1{i} + slask;
	 end
      end
      
      k1 = cell(1,ml);
      for i = 1:ml
	 vik = V{i}(ii,:);
	 slask = 0;
	 for j = 1:mp
	    slask = slask + vik*N^-1*yp(:,j);
	 end
	 k1{i} = slask;
      end
      
      A2 = cell(1,ml);
      for i = 1:ml 
	 A2{i} = zeros(4,4);
	 slask = zeros(4,4);
	 vik = V{i}(ii,:);
	 for j = 1:ml
	    a = (vik*N^-1*V{j}'*V{j})' ;
	    b = xl(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A2{i} = A2{i} + slask;
	 end
      end
      
      k2 = cell(1,ml);
      for i = 1:ml
	 slask = 0;
	 vik = V{i}(ii,:);
	 for j = 1:ml
	    slask = slask + vik*N^-1*V{j}'*V{j}*yl(:,j);
	 end
	 k2{i} = slask;
      end
      
      A3 = cell(1,ml);
      for i = 1:ml
	 vik = V{i}(ii,:);
	 slask = zeros(4,4); 
	 A3{i} = zeros(4,4);
	 for j = 1:mpi
	    a = (vik*N^-1*n(:,j)*n(:,j'))' ;
	    b = xpi(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A3{i} = A3{i} + slask;
	 end
      end
      
      k3 = cell(1,mpi);
      for i = 1:ml
	 slask = 0;
	 vik = V{i}(ii,:);
	 for j = 1:mpi
	    slask = slask + vik*N^-1*n(:,j)*n(:,j)'*ypi(:,j);
	 end
	 k3{i} = slask;
      end
      
      for i = 1:ml
	 vik = V{i}(ii,:);
	 B{end+1} = A{i}-A1{i}-A2{i}-A3{i};
	 k{end+1} = k1{i} + k2{i} + k3{i} - vik*yl(:,i);
      end
   end
end


% Ber�kna m�lfunktion (punkt-punkt delen) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mp >0
   A = cell(1,mp);
   I = eye(3);
   for ii = 1:3
      for i = 1:mp
	 ek = I(ii,:);
	 a = ek';
	 b = xp(:,i);
	 A{i} = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		 a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		 
		 a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		 a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		 
		 a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		 -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		 
		 a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		 a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
      end
      
      
      A1 = cell(1,mp);
      for i = 1:mp
	 ek = I(ii,:);
	 A1{i} = zeros(4,4);
	 slask = zeros(4,4);
	 for j = 1:mp
	    a = (ek*N^-1)' ;
	    b = xp(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A1{i} = A1{i} + slask;
	 end
      end
      
      k1 = cell(1,mp);
      for i = 1:mp
	 ek = I(ii,:);
	 slask = 0;
	 for j = 1:mp
	    slask = slask + ek*N^-1*yp(:,j);
	 end
	 k1{i} = slask;
      end
      
      A2 = cell(1,mp);
      for i = 1:mp 
	 A2{i} = zeros(4,4);
	 slask = zeros(4,4);
	 ek = I(ii,:);
	 for j = 1:ml
	    a = (ek*N^-1*V{j}'*V{j})' ;
	    b = xl(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A2{i} = A2{i} + slask;
	 end
      end
      
      k2 = cell(1,mp);
      for i = 1:mp
	 slask = 0;
	 ek = I(ii,:);
	 for j = 1:ml
	    slask = slask + ek*N^-1*V{j}'*V{j}*yl(:,j);
	 end
	 k2{i} = slask;
      end
      
      A3 = cell(1,mp);
      for i = 1:mp
	 ek = I(ii,:);
	 slask = zeros(4,4); 
	 A3{i} = zeros(4,4);
	 for j = 1:mpi
	    a = (ek*N^-1*n(:,j)*n(:,j'))' ;
	    b = xpi(:,j);
	    slask = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)  a(2)*b(1)+a(1)*b(2) ...
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)-a(2)*b(3); ...
		     
		     a(2)*b(1)+a(1)*b(2)  -a(1)*b(1)+a(2)*b(2)-a(3)*b(3) ...
		     a(3)*b(2)+a(2)*b(3)  -a(3)*b(1)+a(1)*b(3); ...
		     
		     a(3)*b(1)+a(1)*b(3)  a(3)*b(2)+a(2)*b(3) ...
		     -a(1)*b(1)-a(2)*b(2)+a(3)*b(3)  a(2)*b(1)-a(1)*b(2); ...
		     
		     a(3)*b(2)-a(2)*b(3)  -a(3)*b(1)+a(1)*b(3) ...
		     a(2)*b(1)-a(1)*b(2)  a(1)*b(1)+a(2)*b(2)+a(3)*b(3)];
	    A3{i} = A3{i} + slask;
	 end
      end
      
      k3 = cell(1,mp);
      for i = 1:mp
	 slask = 0;
	 ek = I(ii,:);
	 for j = 1:mpi
	    slask = slask + ek*N^-1*n(:,j)*n(:,j)'*ypi(:,j);
	 end
	 k3{i} = slask;
      end
      
      for i = 1:mp
	 ek = I(ii,:);
	 B{end+1} = A{i}-A1{i}-A2{i}-A3{i};
	 k{end+1} = k1{i} + k2{i} + k3{i} - ek*yp(:,i);
      end
   end
end


m = 3*ml+3*mp+mpi;
