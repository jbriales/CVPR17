% This is a script to be called from the outer function BnB_solve

%ordning alpha, q1, q2, q3, q4, s1 = q1^2, s2 = q1q2, s3 = q1q3, s4 = q1q4,
%s5 = q2^2, s6 = q2q3, s7 = q2q4, s8 = q3^2, s9 = q3q4, s10 = q4^2

%konkrav residualer
r = zeros(m,10);
c1 = zeros(m+1,1);
for i = 1:m
   b = B{i};
   ri = [b(1,1) 2*b(1,2) 2*b(1,3) 2*b(1,4) b(2,2) 2*b(2,3) 2*b(2,4) ...
	 b(3,3) 2*b(3,4) b(4,4)];
   r(i,:) = ri;
   c1(i+1) = k{i};
end

At1 = -[1 zeros(1,14); ...
      zeros(m,5) r];

%Affina krav q1-q4 (omr�desgr�nser)
c2 = [-qL(1) qU(1) -qL(2) qU(2) -qL(3) qU(3) -qL(4) qU(4)]';
At2 = zeros(8, 15);
At2(1,2) = -1;
At2(2,2) = 1;
At2(3,3) = -1;
At2(4,3) = 1;
At2(5,4) = -1;
At2(6,4) = 1;
At2(7,5) = -1;
At2(8,5) = 1;

%kon- och affina- krav p� s1-s10
%konvilkor s1 >= q1^2
c3 = [1/4 0 -1/4]';
At3 = zeros(3,15);
At3(1,6) = -1;,
At3(2,2) = -1;
At3(3,6) = -1;

%affint vilkor ist�llet f�r s1 <= q1^2, a*q1+b-s1 >= 0
a = (qU(1)^2-qL(1)^2)/(qU(1) - qL(1));
b = -qL(1)*a + qL(1)^2;
c4 = [b];
At4 = zeros(1,15);
At4(2) = -a;
At4(6) = 1;

%affina villkor ist�llet f�r s2 >= q1q2
c5 = [qU(1)*qU(2) ; qL(1)*qL(2)];
At5 = zeros(2,15);
At5(1,2) = qU(2);
At5(1,3) = qU(1);
At5(1,7) = -1;
At5(2,2) = qL(2);
At5(2,3) = qL(1);
At5(2,7) = -1;

%affina vilkor ist�llet f�r s2 <= q1q2
c6 = [-qU(1)*qL(2); -qL(1)*qU(2)];
At6 = zeros(2,15);
At6(1,2) = -qL(2);
At6(1,3) = -qU(1);
At6(1,7) = 1;
At6(2,2) = -qU(2);
At6(2,3) = -qL(1);
At6(2,7) = 1;

%s3 >= q1q3
c7 = [qU(1)*qU(3) ; qL(1)*qL(3)];
At7 = zeros(2,15);
At7(1,2) = qU(3);
At7(1,4) = qU(1);
At7(1,8) = -1;
At7(2,2) = qL(3);
At7(2,4) = qL(1);
At7(2,8) = -1;

%s3 <= q1q3
c8 = [-qU(1)*qL(3); -qL(1)*qU(3)];
At8 = zeros(2,15);
At8(1,2) = -qL(3);
At8(1,4) = -qU(1);
At8(1,8) = 1;
At8(2,2) = -qU(3);
At8(2,4) = -qL(1);
At8(2,8) = 1;

%s4 >= q1q4
c9 = [qU(1)*qU(4) ; qL(1)*qL(4)];
At9 = zeros(2,15);
At9(1,2) = qU(4);
At9(1,5) = qU(1);
At9(1,9) = -1;
At9(2,2) = qL(4);
At9(2,5) = qL(1);
At9(2,9) = -1;

%s4 <= q1q4
c10 = [-qU(1)*qL(4); -qL(1)*qU(4)];
At10 = zeros(2,15);
At10(1,2) = -qL(4);
At10(1,5) = -qU(1);
At10(1,9) = 1;
At10(2,2) = -qU(4);
At10(2,5) = -qL(1);
At10(2,9) = 1;

%konvilkor s5 >= q2^2
c11 = [1/4 0 -1/4]';
At11 = zeros(3,15);
At11(1,10) = -1;,
At11(2,3) = -1;
At11(3,10) = -1;

%affint vilkor a*q2+b-s5 >= 0 ist�llet f�r s5 <= q2^2.
a = (qU(2)^2-qL(2)^2)/(qU(2) - qL(2));
b = -qL(2)*a + qL(2)^2;
c12 = [b];
At12 = zeros(1,15);
At12(3) = -a;
At12(10) = 1;

%s6 >= q2q3
c13 = [qU(2)*qU(3) ; qL(2)*qL(3)];
At13 = zeros(2,15);
At13(1,3) = qU(3);
At13(1,4) = qU(2);
At13(1,11) = -1;
At13(2,3) = qL(3);
At13(2,4) = qL(2);
At13(2,11) = -1;

%s6 <= q2q3
c14 = [-qU(2)*qL(3); -qL(2)*qU(3)];
At14 = zeros(2,15);
At14(1,3) = -qL(3);
At14(1,4) = -qU(2);
At14(1,11) = 1;
At14(2,3) = -qU(3);
At14(2,4) = -qL(2);
At14(2,11) = 1;

%s7 >= q2q4
c15 = [qU(2)*qU(4) ; qL(2)*qL(4)];
At15 = zeros(2,15);
At15(1,3) = qU(4);
At15(1,5) = qU(2);
At15(1,12) = -1;
At15(2,3) = qL(4);
At15(2,5) = qL(2);
At15(2,12) = -1;

%s7 <= q2q4
c16 = [-qU(2)*qL(4); -qL(2)*qU(4)];
At16 = zeros(2,15);
At16(1,3) = -qL(4);
At16(1,5) = -qU(2);
At16(1,12) = 1;
At16(2,3) = -qU(4);
At16(2,5) = -qL(2);
At16(2,12) = 1;

%konvilkor s8 >= q3^2
c17 = [1/4 0 -1/4]';
At17 = zeros(3,15);
At17(1,13) = -1;,
At17(2,4) = -1;
At17(3,13) = -1;

%affint vilkor a*q3+b-s8 >= 0 ist�llet f�r s8 <= q3^2
a = (qU(3)^2-qL(3)^2)/(qU(3) - qL(3));
b = -qL(3)*a + qL(3)^2;
c18 = [b];
At18 = zeros(1,15);
At18(4) = -a;
At18(13) = 1;

%s9 >= q3q4
c19 = [qU(3)*qU(4) ; qL(3)*qL(4)];
At19 = zeros(2,15);
At19(1,4) = qU(4);
At19(1,5) = qU(3);
At19(1,14) = -1;
At19(2,4) = qL(4);
At19(2,5) = qL(3);
At19(2,14) = -1;

%s9 <= q3q4
c20 = [-qU(3)*qL(4); -qL(3)*qU(4)];
At20 = zeros(2,15);
At20(1,4) = -qL(4);
At20(1,5) = -qU(3);
At20(1,14) = 1;
At20(2,4) = -qU(4);
At20(2,5) = -qL(3);
At20(2,14) = 1;

%konvilkor s10 >= q4^2
c21 = [1/4 0 -1/4]';
At21 = zeros(3,15);
At21(1,15) = -1;,
At21(2,5) = -1;
At21(3,15) = -1;

%affint vilkor a*q4+b-s10 >= 0 ist�llet f�r s10 <= q4^2
a = (qU(4)^2-qL(4)^2)/(qU(4) - qL(4));
b = -qL(4)*a + qL(4)^2;
c22 = [b];
At22 = zeros(1,15);
At22(5) = -a;
At22(15) = 1;

%Euclidiska villkor
At23 = zeros(2,15);
At23(1,6) = -1;
At23(1,10) = -1;
At23(1,13) = -1;
At23(1,15) = -1;
At23(2,6) = 1;
At23(2,10) = 1;
At23(2,13) = 1;
At23(2,15) = 1;
c23 = [-1 1]';

%Samla affina vilkor
At_affin = [At2; At4; At5; At6; At7; At8; At9; At10; At12; At13; At14; At15; ...
      At16; At18; At19; At20; At22; At23];

c_affin = [c2; c4; c5; c6; c7; c8; c9; c10; c12; c13; c14; c15; ...
      c16; c18; c19; c20; c22; c23];

%m�lfunktion -alpha
b = [-1; zeros(14,1)];


%Samla alla vilkor
At = [At_affin; At1; At3; At11; At17; At21];
c = [c_affin; c1; c3; c11; c17; c21];
K.l = size(At_affin, 1);
K.q = [size(At1,1)  size(At3,1) size(At11,1) size(At17,1) size(At21,1)];