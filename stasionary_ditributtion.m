clear;
clc;
D0 = [-1 0.5; 1 -3];
D1 = [0 0.5; 1 1];
I = [1 0; 0 1];
E = [1; 1];
E3 = [1; 1; 1];
miu = 2;
A0 = D1;
A2 = miu*I;
A1 = D0 - A2;
R0 = [0 0; 0 0];
Y = inv(A1);
R1 = -A0*Y;
counter = 1;
while norm(R1-R0) ~= 0
    R0 = R1;
    R1 = -(A0 + R0*R0*A2)*Y;
    counter = counter + 1 ;
end
IN = inv(I-R1);
A = D0 + R1*A2;
B = IN*E; 
p0 =[0; 0];
C(1,1) = A(1,1);
C(1,2) = A(2,1);
C(2,1) = B(1,1);
C(2,2) = B(2,1);
p0 = linsolve(C,[0;1]);
p0 = transpose(p0);

 % form teacher;
lamda_out = (p0*(R1*IN*A2))*E; 

% this formulas form David Anthony Green
x0 = inv(p0*IN*A0*E)*p0*R1*A2;
% set the Q0
Q0(3,3) = 0;
for i = 1:2
    for j = 1:2
        Q0(i,j) = D0(i,j);
    end
end
QQ0 = D1*E;
for i = 1:2
    Q0(i,3) = QQ0(i);
end
Q0(3,3) = -miu;
% set the Q1
Q1(3,3) = 0;
QQ1 = miu*x0;
for j = 1:2
    Q1(3,j) = QQ1(j);
end
Q1(3,3) = (1-x0*E)*miu;
% compute the static of output
P = inv(-Q0)*Q1;
P1 = transpose(P);
M(4,3) = 0;
for i = 1:3
    for j = 1:3
        M(i,j) = P1(i,j);
    end
end
for i = 1:3
    M(i,i) = M(i,i)-1;
end
for i = 1:3
    M(4,i) = 1;
end
Pi = linsolve(M,[0; 0; 0; 1]);
Pi = transpose(Pi);
lamda_o =1/(Pi*inv(-Q0)*E3);

Ex = 1/lamda_o;
EEx = Ex*Ex;
Exx = 2*Pi*inv(-Q0)*inv(-Q0)*E3;
variance = Exx - EEx ;








