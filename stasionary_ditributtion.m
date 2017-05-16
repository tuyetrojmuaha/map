clear;
clc;
D0 = [-0.1 0.05; 0.1 -0.3];
D1 = [0 0.05; 0.1 0.1];
k = 1;
D0 = k.*D0;
D1 = k.*D1;
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
%p0 =[0; 0];
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
stddev = sqrt(variance);

alpha = lamda_o*Pi*inv(-Q0);
variance1 = 2*(1/lamda_o)*alpha*inv(-Q0)*E3 - 1/(lamda_o*lamda_o);

% compute the correlations;

t = 0;
a = lamda_o^2*Pi*(inv(-Q0))*P^1*(inv(-Q0))*E3 - 1 ;
b = 2*lamda_o^2*Pi*(inv(-Q0))*(inv(-Q0))*E3 - 1 ;
correlations_of_t =  a/b ;

% compute the approximations of MAP(l);

l = 10;
Ql0 = zeros(l*2,l*2);
Ql1 = zeros(l*2,l*2);

for i = 1:2
    for j = 1:2
        Ql0(i,j) = D0(i,j);
    end
end
for m = 1 : (l-1)
    for i = 1:2
        for j = 1:2
        Ql0((m-1)*2+i,m*2+j) = A0(i,j);
        end
    end
    for i = 1:2
        for j = 1:2
        Ql0(m*2+i,m*2+j) = A1(i,j);
        end
    end
end
for m = 1 : (l-1)
    for i = 1:2
        for j = 1:2
        Ql1(m*2+i,(m-1)*2+j) = A2(i,j);
        end
    end
end

Il = zeros(l*2,l*2);
for i = 1:l*2
    Il(i,i) = 1;
end
El = ones(l*2,1);
Pl = (inv(-Ql0))*Ql1;
Pl1 = Pl - Il;
Pl2 = transpose(Pl1);
Ml = ones(l*2+1,l*2);
for i = 1:2*l
    for j = 1:2*l
        Ml(i,j) = Pl1(i,j);
    end
end
ll = ones(2*l+1,1);
ll(2*m+1,1) = 1;
pil = linsolve(Ml,ll);
pil = transpose(pil);
lamdal = 1/(pil*(inv(-Ql0))*El);
tl = 0;
al = lamdal^2*pil*(inv(-Ql0))*Pl^tl*(inv(-Ql0))*El - 1 ;
bl = 2*lamdal^2*pil*(inv(-Ql0))*(inv(-Ql0))*El - 1 ;
correlations_of_tl =  al/bl ;

















