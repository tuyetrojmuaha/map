clear;
clc;
D0 = [-0.1 0.05; 0.1 -0.3];
D1 = [0 0.05; 0.1 0.1];
k = 95;
D0 = k*D0;
D1 = k*D1;
I = [1 0; 0 1];
E = [1; 1];
miu = 10;
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













