clear;
clc;
D0 = [-0.1 0.05; 0.1 -0.3];
D1 = [0 0.05; 0.1 0.1];
D1 = 100*D1;
D0 = 100*D0;
I = [1 0; 0 1];
E = [1; 1];
miu = 2;
epsilon = 0.0000000001;
A0 = D1;
A2 = miu.*I;
A1 = D0 - A2;
R0 = [0 0; 0 0];
Y = inv(A1);
R1 = -A0.*Y;
counter = 1;
while norm(R1-R0) ~= 0
    R0 = R1;
    R1 = -(A0 + R0.*R0.*A2).*Y;
    counter = counter + 1 ;
end
% for i=1:2
%     R0 = R1;
%     R1 = -(A0 + R0.*R0.*A2).*Y;
% end

A = inv(I-R1)*E;
B = D0 + R1.*A2;









