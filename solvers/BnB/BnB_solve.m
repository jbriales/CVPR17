function [Rmin,tmin] = BnB_solve(correspondences,tol)
% [Rmin,tmin] = BnB_solve(B,k,m)
% Solve the optimal R and t from the data given by BnB_setData

% take correspondences and create variables used in BnB solver
% BnB_setData(correspondences)
run BnB_setData

% optimset fmincon;
% options = ans;
% options = optimset(options,'TolFun', 1e-16,'TolCon', 1e-16,'Display','off','LargeScale','off','DiffMinChange',1e-16);

global qL;
global qU;
%global Bconv;
%global cconv;
%global kconv;
% $$$ global Bconv1;
% $$$ global cconv1;
% $$$ global kconv1;
% $$$ global Bconv2;
% $$$ global cconv2;
% $$$ global kconv2;

qL = [0 -1 -1 -1]';
qU = [1 1 1 1]';

qLU = [qL qU];
newinterv = cell(1,1);
newinterv{1} = qLU;
%%%%%%%%%%%%%%%%%%%
intervnr = [];
mf = [];
lokmin = [];
ints = [];
%%%%%%%%%%%%%%%%%%

if nargin < 2
  tol = 1e-10;
  warning('No tolerance given, default used: %f',tol)
end
  
minf = inf;
intervsize = 8;
%while (minf > tol) || (intervsize > tol);
while (intervsize > tol);
  interv = newinterv;
  %%%%%%%%%%%%%%%%%%%%%%
  intervnr = [intervnr; length(interv)];
  ints = [ints; intervsize];
  %%%%%%%%%%%%%%%%%%%%%%
  newinterv = cell(0);
  disp(strcat('Number range:',num2str(length(interv))));
  disp(strcat('Interval size:',num2str(intervsize)));
  disp(strcat('Minimum value:',num2str(minf)));
  qL = interv{1}(:,1);
  qU = interv{1}(:,2);
  for i = 1:length(interv)
    qL = interv{i}(:,1);
    qU = interv{i}(:,2);
    
    BnB_genSeDuMiProbEucl
    clear pars;
    pars.fid=0;
    pars.eps = 0;
    [x,y,info] = sedumi(At,b,c,K,pars);
    
    qconv = y(2:5);
    
    fconv = y(1)^2;
    if info.dinf == 1
      disp('infeasible');
      fconv = inf;
    end
    
    if qconv ~= 0
      qconv = qconv./norm(qconv);
    else
      qconv = [1 1 1 1]'./norm([1 1 1 1]');
    end
    
    f = BnB_goalfun(qconv);
    
    % $$$      if norm(qU-qL) > 0.5;
    % $$$ 	 [qlok,flok] = fmincon(@goalfun,qconv,[],[],[],[],qL,qU,[], ...
    % $$$ 			      options);
    % $$$ 	 if flok < minf;
    % $$$ 	    minf = flok;
    % $$$ 	    qmin = qlok;
    % $$$ 	 end
    % $$$      end
    
    if f < minf;
      minf = f;
      qmin = qconv;
    end
    
    if fconv <= minf;
      slask = qU-qL;
      ind = find(max(slask) == slask);
      ind(1);
      qL1 = qL;
      qU1 = qU;
      qL2 = qL;
      qU2 = qU;
      qU1(ind(1)) = (qU(ind(1))+qL(ind(1)))/2;
      qL2(ind(1)) = (qU(ind(1))+qL(ind(1)))/2;
      newinterv{end+1} = [qL1 qU1];
      newinterv{end+1} = [qL2 qU2];
    end
  end
  if length(newinterv) > 0;
    intervsize = prod((newinterv{1}(:,1) - newinterv{1}(:,2)))* ...
      length(newinterv);
  else
    intervsize = 0;
  end
end

xx = qmin(1);
yy = qmin(2);
zz = qmin(3);
ww = qmin(4);

Rmin = [...
  xx^2-yy^2-zz^2+ww^2 2*xx*yy-2*ww*zz 2*zz*xx+2*ww*yy; ...
  2*xx*yy+2*ww*zz -xx^2+yy^2-zz^2+ww^2 2*yy*zz-2*ww*xx; ...
  2*zz*xx-2*ww*yy 2*yy*zz+2*ww*xx -xx^2-yy^2+zz^2+ww^2];

tmin = 0;
for i = 1:mp
  tmin = tmin - Rmin*xp(:,i) + yp(:,i);
end
for i = 1:ml
  tmin = tmin - V{i}'*V{i}*Rmin*xl(:,i) + V{i}'*V{i}*yl(:,i);
end
for i= 1:mpi
  tmin = tmin - n(:,i)*n(:,i)'*Rmin*xpi(:,i) + n(:,i)*n(:,i)'*ypi(:,i);
end
tmin = N^-1*tmin;

end