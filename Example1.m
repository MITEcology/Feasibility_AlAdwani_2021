clear all; close all; clc;

%% Part 1

syms r1 a11 a12 b1 ...
     r2 a21 a22 b2 ...
     N1 N2 f1 f2 ...
     x y 
     
F1=@(n1,n2) r1+a11*n1+a12*n2+b1*n1*n2;
F2=@(n1,n2) r2+a21*n1+a22*n2+b2*n1*n2;

ResN1=@(N1,N2) det([a11 + b1*N2, f1; ...
                    a21 + b2*N2, f2]);
     
ResN2=@(N1,N2) det([a12 + b1*N1, f1; ...
                    a22 + b2*N1, f2]);
 
T=@(N1,N2) det([a21 + b2*N2, a11 + b1*N2; ...
                a22 + b2*N1, a12 + b1*N1]);
 
X=@(N1,N2) det([a11 + b1*N2, F1(N1,N2); ...
                a21 + b2*N2, F2(N1,N2)]);
     
Y=@(N1,N2) det([a12 + b1*N1, F1(N1,N2); ...
                a22 + b2*N1, F2(N1,N2)]);
     
J=@(N1,N2) det([a11 + b1*N2, a12 + b1*N1; ...
                a21 + b2*N2, a22 + b2*N1]);
     
G=@(N1,N2) simplify(T(N1,N2)*J(N1,N2)/(X(N1,N2)*Y(N1,N2))); 

H = G(1/x,1/y);

LL = simplify(taylor(H, [x,y], [0,0],'Order', 9));

Func = matlabFunction(LL);

Sr =@(x,y) Func(a11,a12,a21,a22,b1,b2,r1,r2,x,y);
 
[cx,tx] = coeffs(Sr,[x,y]);

Sigma_00 = cx(tx == x*y);
[Sigma_00U,Sigma_00D] = numden(Sigma_00);

Sigma_10 = cx(tx == x^2*y);
[Sigma_10U,Sigma_10D] = numden(Sigma_10);

Sigma_01 = cx(tx == x*y^2);
[Sigma_01U,Sigma_01D] = numden(Sigma_01);

Sigma_11 = cx(tx == x^2*y^2);
[Sigma_11U,Sigma_11D] = numden(Sigma_11);

Sigma_20 = cx(tx == x^3*y);
[Sigma_20U,Sigma_20D] = numden(Sigma_20);

Sigma_02 = cx(tx == x*y^3);
[Sigma_02U,Sigma_02D] = numden(Sigma_02);

Sigma_21 = cx(tx == x^3*y^2);
[Sigma_21U,Sigma_21D] = numden(Sigma_21);

Sigma_12 = cx(tx == x^2*y^3);
[Sigma_12U,Sigma_12D] = numden(Sigma_12);

Sigma_30 = cx(tx == x^4*y);
[Sigma_30U,Sigma_30D] = numden(Sigma_30);

Sigma_03 = cx(tx == x*y^4);
[Sigma_03U,Sigma_03D] = numden(Sigma_03);

Sigma_31 = cx(tx == x^4*y^2);
[Sigma_31U,Sigma_31D] = numden(Sigma_31);

Sigma_13 = cx(tx == x^2*y^4);
[Sigma_13U,Sigma_13D] = numden(Sigma_13);

Sigma_22 = cx(tx == x^3*y^3);
[Sigma_22U,Sigma_22D] = numden(Sigma_22);

%syms x1 y1 x2 y2 x1 c1 c2 s1 s2 
%S = collect([1, 1; c1*x1+c2*y1, c1*x2+c2*y2]*[(x1-s1)*(y1-s2),0;0,(x2-s1)*(y2-s2)]* ...
%            [1, c1*x1+c2*y1; 1, c1*x2+c2*y2], [s1,s2, c1,c2]);

%% Part 2

syms s1 s2 ...
     sigma_10 sigma_01 sigma_11 sigma_20 sigma_21 sigma_02 ...
     sigma_12 sigma_22 sigma_30 sigma_31 sigma_03 sigma_13 ...
     c1 c2 

S_sigma = @(s1,s2) [2*s1*s2 - sigma_01*s1 - sigma_10*s2 + sigma_11, ...
                   c1*(sigma_10*s1*s2 - sigma_11*s1 - sigma_20*s2 + sigma_21) + ...
                   c2*(sigma_01*s1*s2 - sigma_11*s2 - sigma_02*s1 + sigma_12); ...
                   c1*(sigma_10*s1*s2 - sigma_11*s1 - sigma_20*s2 + sigma_21) + ...
                   c2*(sigma_01*s1*s2 - sigma_11*s2 - sigma_02*s1 + sigma_12), ...
                   c1^2*(sigma_20*s1*s2 - sigma_21*s1 - sigma_30*s2 + sigma_31) + ...
                   2*c1*c2*(sigma_11*s1*s2 - sigma_12*s1 - sigma_21*s2 + sigma_22) + ...
                   c2^2*(sigma_02*s1*s2 - sigma_12*s2 - sigma_03*s1 + sigma_13)]; 
        
Trace_S = @(s1,s2) trace(S_sigma(s1,s2)); 
det_S = @(s1,s2) det(S_sigma(s1,s2)); 
 
[cT00,tT00] = coeffs(Trace_S(s1,s2),[s1,s2]);
[cD00,tD00] = coeffs(det_S(s1,s2),[s1,s2]);
mtr00 = -sign(cT00(tT00 == 1));
det00 =  sign(cD00(tD00 == 1)); 

[Ti0_N,Ti0_D] = numden(Trace_S(1/s1,0));
[Di0_N,Di0_D] = numden(det_S(1/s1,0));
[cTi0,tTi0] = coeffs(Ti0_N,s1);
[cDi0,tDi0] = coeffs(Di0_N,s1);
mtri0 = -sign(cTi0(tTi0 == 1));
deti0 =  sign(cDi0(tDi0 == 1)); 

[T0i_N,T0i_D] = numden(Trace_S(0,1/s2));
[D0i_N,D0i_D] = numden(det_S(0,1/s2));
[cT0i,tT0i] = coeffs(T0i_N,s2);
[cD0i,tD0i] = coeffs(D0i_N,s2);
mtr0i = -sign(cT0i(tT0i == 1));
det0i =  sign(cD0i(tD0i == 1)); 

[Tii_N,Tii_D] = numden(Trace_S(1/s1,1/s2));
[Dii_N,Dii_D] = numden(det_S(1/s1,1/s2));
[cTii,tTii] = coeffs(Tii_N,[s1,s2]);
[cDii,tDii] = coeffs(Dii_N,[s1,s2]);
mtrii = -sign(cTii(tTii == 1));
detii =  sign(cDii(tDii == 1)); 

V00 = (1 - mtr00)/2 + (1 - mtr00 * det00)/2;
Vi0 = (1 - mtri0)/2 + (1 - mtri0 * deti0)/2;
V0i = (1 - mtr0i)/2 + (1 - mtr0i * det0i)/2;
Vii = (1 - mtrii)/2 + (1 - mtrii * detii)/2;

NpRoots = (V00 - Vi0 + Vii - V0i)/2; 

%% Part 3

sigma_10 = Sigma_10;  sigma_01 = Sigma_01;  sigma_11 = Sigma_11;
sigma_20 = Sigma_20;  sigma_21 = Sigma_21;  sigma_02 = Sigma_02;  
sigma_12 = Sigma_12;  sigma_22 = Sigma_22;  sigma_30 = Sigma_30; 
sigma_31 = Sigma_31;  sigma_03 = Sigma_03;  sigma_13 = Sigma_13;

F_mtr00 = matlabFunction(subs(mtr00));
F_mtri0 = matlabFunction(subs(mtri0));
F_mtr0i = matlabFunction(subs(mtr0i));
F_mtrii = matlabFunction(subs(mtrii));

F_det00 = matlabFunction(subs(det00));
F_deti0 = matlabFunction(subs(deti0));
F_det0i = matlabFunction(subs(det0i));
F_detii = matlabFunction(subs(detii));

F_NpRoots = matlabFunction(subs(NpRoots));

%% Part 4

clc; 

sample = 10^5;

a11 = 20*(rand(sample,1)-0.5); 
a12 = 20*(rand(sample,1)-0.5);
a21 = 20*(rand(sample,1)-0.5); 
a22 = 20*(rand(sample,1)-0.5); 
b1  = 20*(rand(sample,1)-0.5); 
b2  = 20*(rand(sample,1)-0.5);
c1 = 1*ones(sample,1); 
c2 = 0*ones(sample,1);
r1  = 20*(rand(sample,1)-0.5); 
r2  = 20*(rand(sample,1)-0.5);

Signs = [F_mtr00(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_det00(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_mtri0(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_deti0(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_mtr0i(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_det0i(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_mtrii(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_detii(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2), ...
         F_NpRoots(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2)]; 

unique_Signs = unique(Signs,'rows')
[rows,~] = find(unique_Signs(:,end)>0);
unique_nonzero_Signs = sortrows(unique_Signs(rows,1:end-1))

% Part 5: Checking independency on c1 and c2

cc1 = rand*ones(sample,1); 
cc2 = rand*ones(sample,1);

Count1 = F_NpRoots(a11,a12,a21,a22,b1,b2,c1,c2,r1,r2); 
Count2 = F_NpRoots(a11,a12,a21,a22,b1,b2,cc1,cc2,r1,r2); 
K = Count1 - Count2; 
KK = find(abs(K)>0);
check = length(KK)

% Calculating Probability of feasibility 

probability = sum(logical(Count1))/(sample)