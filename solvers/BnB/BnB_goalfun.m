function y = BnB_goalfun(q);

global B;
global k;
global m;

y = 0;
for i = 1:m;
   y = y + (q'*B{i}*q + k{i})^2;
end
