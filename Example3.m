clear all; close all; clc; 

%% Part 1

syms r1 r2 r3 ...
     N1 N2 N3 ...
     F1 F2 F3 ... 
     a11 a12 a13 b1 ... 
     a21 a22 a23 b2 ...
     a31 a32 a33 b3 
 
f1 = @(N1,N2,N3) r1+a11*N1+a12*N2+a13*N3+b1*N1*N2;
f2 = @(N1,N2,N3) r2+a21*N1+a22*N2+a23*N3+b2*N1*N2;
f3 = @(N1,N2,N3) r3+a31*N1+a32*N2+a33*N3+b3*N1*N2;

% Step 1a: Obtaining Res_N2N3(N1)

l1 = 3; 

H1row = [f1(N1,N2,N3),  f2(N1,N2,N3), f3(N1,N2,N3)]; 
H1col = [1, N2, N3]; 

for i = 1:l1
    
    [c,T] = coeffs(H1row(i),[N2,N3]);
    
    for j = 1:l1
        
        if isempty(c(T==H1col(j))) 
            
            M_N1(i,j) = 0;
            
        else
            
            M_N1(i,j) = c(T==H1col(j)); 
            
        end
        
    end
    
end

Res_N2N3 = collect(det(M_N1),N1);
[cRes_N2N3,tRes_N2N3] = coeffs(Res_N2N3,N1); 

Mp_N1 = M_N1;
Mp_N1(:,1) = [F1, F2, F3]; 

T1 = collect(det(Mp_N1),[F1 F2 F3]);
T11 = subs(T1, [F1 F2 F3], [1 0 0]); 
T12 = subs(T1, [F1 F2 F3], [0 1 0]); 
T13 = subs(T1, [F1 F2 F3], [0 0 1]); 

t11 = matlabFunction(T11); 
fT11 =@(N1)    t11(N1,a22,a23,a32,a33,b2,b3);
t12 = matlabFunction(T12); 
fT12 =@(N1)    t12(N1,a12,a13,a32,a33,b1,b3);
t13 = matlabFunction(T13); 
fT13 =@(N1)    t13(N1,a12,a13,a22,a23,b1,b2);

% Step 1b: Obtaining Res_N1N3(N2)

l2 = 3; 

H2row = [f1(N1,N2,N3),  f2(N1,N2,N3), f3(N1,N2,N3)]; 
H2col = [1, N1, N3];

for i = 1:l2
    
    [c,T] = coeffs(H2row(i),[N1,N3]);
    
    for j = 1:l2
        
        if isempty(c(T==H2col(j))) 
            
            M_N2(i,j) = 0;
            
        else
            
            M_N2(i,j) = c(T==H2col(j)); 
            
        end
        
    end
    
end

Res_N1N3 = collect(det(M_N2),N2);
[cRes_N1N3,tRes_N1N3] = coeffs(Res_N1N3,N2); 

Mp_N2 = M_N2;
Mp_N2(:,1) = [F1, F2, F3]; 

T2 = collect(det(Mp_N2),[F1 F2 F3]);
T21 = subs(T2, [F1 F2 F3], [1 0 0]); 
T22 = subs(T2, [F1 F2 F3], [0 1 0]); 
T23 = subs(T2, [F1 F2 F3], [0 0 1]); 

t21 = matlabFunction(T21); 
fT21 =@(N2)    t21(N2,a21,a23,a31,a33,b2,b3);
t22 = matlabFunction(T22); 
fT22 =@(N2)    t22(N2,a11,a13,a31,a33,b1,b3);
t23 = matlabFunction(T23); 
fT23 =@(N2)    t23(N2,a11,a13,a21,a23,b1,b2);

% Step 1c: Obtaining Res_N1N2(N3)

l3 = 15; 

H3row = [f1(N1,N2,N3) ,N1*f1(N1,N2,N3) ,N2*f1(N1,N2,N3) ,N1*N2*f1(N1,N2,N3) ,N1^2*f1(N1,N2,N3) ,N2^2*f1(N1,N2,N3),...
         N1*f2(N1,N2,N3) ,N2*f2(N1,N2,N3) ,N1*N2*f2(N1,N2,N3) ,N1^2*f2(N1,N2,N3) ,N2^2*f2(N1,N2,N3) ,N1*f3(N1,N2,N3),...
         N2*f3(N1,N2,N3), N1*N2*f3(N1,N2,N3),N2^2*f3(N1,N2,N3)]; 
H3col = [1 ,N1 ,N2 ,N1*N2 ,N1^2 ,N2^2 ,N1^3 ,N1^2*N2 ,N1^3*N2 ,N1^4 ,N1^2*N2^2 ,N1*N2^2 ,N2^3 ,N1*N2^3 ,N2^4];

for i = 1:l3
    
    [c,T] = coeffs(H3row(i),[N1,N2]);
    
    for j = 1:l3
        
        if isempty(c(T==H3col(j))) 
            
            M_N3(i,j) = 0;
            
        else
            
            M_N3(i,j) = c(T==H3col(j)); 
            
        end
        
    end
    
end

syms ep; 
Mp3 = M_N3 - ep*eye(l3);
[cMp3,tMp3] = coeffs(det(Mp3),ep);  
Res_N1N2 = cMp3(tMp3==ep^3); 
[cRes_N1N2,tRes_N1N2] = coeffs(Res_N1N2,N3); 

Mp_N3 = M_N3 - ep*eye(l3);
Mp_N3(:,1) = [F1 ,N1*F1 ,N2*F1 ,N1*N2*F1 ,N1^2*F1 ,N2^2*F1,...
              N1*F2 ,N2*F2 ,N1*N2*F2 ,N1^2*F2, N2^2*F2 ,N1*F3,...
              N2*F3, N1*N2*F3,N2^2*F3];
                
[cMp_N3,tMp_N3] = coeffs(det(Mp_N3),ep); 
T3 = collect(cMp_N3(tMp_N3==ep^3),[F1 F2 F3]); 
T31 = subs(T3, [F1 F2 F3], [1 0 0]); 
T32 = subs(T3, [F1 F2 F3], [0 1 0]); 
T33 = subs(T3, [F1 F2 F3], [0 0 1]); 

t31 = matlabFunction(T31); 
fT31 =@(N1,N2,N3)    t31(N1,N2,N3,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);
t32 = matlabFunction(T32); 
fT32 =@(N1,N2,N3)    t32(N1,N2,N3,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);
t33 = matlabFunction(T33); 
fT33 =@(N1,N2,N3)    t33(N1,N2,N3,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);


T =@(N1,N2,N3) det([fT11(N1)       fT12(N1)       fT13(N1); ...
                    fT21(N2)       fT22(N2)       fT23(N2); ...
                    fT31(N1,N2,N3) fT32(N1,N2,N3) fT33(N1,N2,N3)]); 
                
Jacob = det(jacobian([f1(N1,N2,N3), f2(N1,N2,N3), f3(N1,N2,N3)], [N1, N2, N3])); 
fJ = matlabFunction(Jacob);
J =@(N1,N2) fJ(N1,N2,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3);

R1 = matlabFunction(Res_N2N3);
Res1 =@(N1) R1(N1,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);
R2 = matlabFunction(Res_N1N3);
Res2 =@(N2) R2(N2,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);
R3 = matlabFunction(Res_N1N2);
Res3 =@(N3) R3(N3,a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);

syms x y z

H1 = 1/Res1(1/x);
L1 = taylor(H1,'Order',9);
LL1 = matlabFunction(L1); 
Resultant1 =@(x) LL1(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3,x);

H2 = 1/Res2(1/y);
L2 = taylor(H2,'Order', 9);
LL2 = matlabFunction(L2); 
Resultant2 =@(y) LL2(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3,y);

H3 = 1/Res3(1/z);
L3 = taylor(H3,'Order', 9);
LL3 = matlabFunction(L3); 
Resultant3 =@(z) LL3(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3,z);

[cTJ, tTJ] = coeffs(T(x,y,z)*J(x,y),[x y z]);
[cx, tx] = coeffs(Resultant1(x),x);
[cy, ty] = coeffs(Resultant2(y),y);
[cz, tz] = coeffs(Resultant3(z),z);

Sigma_000 = 0;
Sigma_001 = 0; Sigma_010 = 0; Sigma_011 = 0; 
Sigma_100 = 0; Sigma_101 = 0; Sigma_110 = 0; 
Sigma_111 = 0; Sigma_200 = 0; Sigma_201 = 0; 
Sigma_210 = 0; Sigma_211 = 0; Sigma_300 = 0; 
Sigma_301 = 0; Sigma_310 = 0; Sigma_311 = 0; 
                                           

for i = 1:length(tTJ)
    
    
    powx = polynomialDegree(tTJ(i),x);
    powy = polynomialDegree(tTJ(i),y);
    powz = polynomialDegree(tTJ(i),z);
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_000 = Sigma_000 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_100 = Sigma_100 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_200 = Sigma_200 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_300 = Sigma_300 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+0)); 
    end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_010 = Sigma_010 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_110 = Sigma_110 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_210 = Sigma_210 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_310 = Sigma_310 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+0)); 
    end
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_001 = Sigma_001 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_101 = Sigma_101 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_201 = Sigma_201 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_301 = Sigma_301 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+1)); 
    end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_011 = Sigma_011 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_111 = Sigma_111 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_211 = Sigma_211 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_311 = Sigma_311 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+1)); 
    end
    
end

%% part 2

syms N11 N12 N21 N22 N31 N32 s1 s2 s3

m=@(N1,N2,N3) [1; N1];

W=[m(N11,N21,N31),m(N12,N22,N32)]; 
D = diag([(N11-s1)*(N21-s2)*(N31-s3),(N12-s1)*(N22-s2)*(N32-s3)]);
S=collect(W*D*transpose(W),[s1 s2 s3]); 

syms  S000 S100 S200 S300 S010 S110 S210 S310 S020 S120 S220 S320 S030 S130 S230 S330 ...
      S001 S101 S201 S301 S011 S111 S211 S311 S021 S121 S221 S321 S031 S131 S231 S331 ...
      S002 S102 S202 S302 S012 S112 S212 S312 S022 S122 S222 S322 S032 S132 S232 S332 ...
      S003 S103 S203 S303 S013 S113 S213 S313 S023 S123 S223 S323 S033 S133 S233 S333
  
sums =[N11^1*N21^0*N31^0 + N12^1*N22^0*N32^0, ...
       N11^2*N21^0*N31^0 + N12^2*N22^0*N32^0, ...
       N11^3*N21^0*N31^0 + N12^3*N22^0*N32^0, ...
       N11^0*N21^1*N31^0 + N12^0*N22^1*N32^0, ...
       N11^1*N21^1*N31^0 + N12^1*N22^1*N32^0, ...
       N11^2*N21^1*N31^0 + N12^2*N22^1*N32^0, ...
       N11^3*N21^1*N31^0 + N12^3*N22^1*N32^0, ...
       N11^0*N21^2*N31^0 + N12^0*N22^2*N32^0, ...
       N11^1*N21^2*N31^0 + N12^1*N22^2*N32^0, ...
       N11^2*N21^2*N31^0 + N12^2*N22^2*N32^0, ...
       N11^3*N21^2*N31^0 + N12^3*N22^2*N32^0, ...
       N11^0*N21^3*N31^0 + N12^0*N22^3*N32^0, ...  
       N11^1*N21^3*N31^0 + N12^1*N22^3*N32^0, ...
       N11^2*N21^3*N31^0 + N12^2*N22^3*N32^0, ...
       N11^3*N21^3*N31^0 + N12^3*N22^3*N32^0, ...   
       N11^0*N21^0*N31^1 + N12^0*N22^0*N32^1, ...
       N11^1*N21^0*N31^1 + N12^1*N22^0*N32^1, ...
       N11^2*N21^0*N31^1 + N12^2*N22^0*N32^1, ...
       N11^3*N21^0*N31^1 + N12^3*N22^0*N32^1, ...
       N11^0*N21^1*N31^1 + N12^0*N22^1*N32^1, ...
       N11^1*N21^1*N31^1 + N12^1*N22^1*N32^1, ...
       N11^2*N21^1*N31^1 + N12^2*N22^1*N32^1, ...
       N11^3*N21^1*N31^1 + N12^3*N22^1*N32^1, ...
       N11^0*N21^2*N31^1 + N12^0*N22^2*N32^1, ...
       N11^1*N21^2*N31^1 + N12^1*N22^2*N32^1, ...
       N11^2*N21^2*N31^1 + N12^2*N22^2*N32^1, ...
       N11^3*N21^2*N31^1 + N12^3*N22^2*N32^1, ...
       N11^0*N21^3*N31^1 + N12^0*N22^3*N32^1, ...  
       N11^1*N21^3*N31^1 + N12^1*N22^3*N32^1, ...
       N11^2*N21^3*N31^1 + N12^2*N22^3*N32^1, ...
       N11^3*N21^3*N31^1 + N12^3*N22^3*N32^1, ...
       N11^0*N21^0*N31^2 + N12^0*N22^0*N32^2, ...
       N11^1*N21^0*N31^2 + N12^1*N22^0*N32^2, ...
       N11^2*N21^0*N31^2 + N12^2*N22^0*N32^2, ...
       N11^3*N21^0*N31^2 + N12^3*N22^0*N32^2, ...
       N11^0*N21^1*N31^2 + N12^0*N22^1*N32^2, ...
       N11^1*N21^1*N31^2 + N12^1*N22^1*N32^2, ...
       N11^2*N21^1*N31^2 + N12^2*N22^1*N32^2, ...
       N11^3*N21^1*N31^2 + N12^3*N22^1*N32^2, ...
       N11^0*N21^2*N31^2 + N12^0*N22^2*N32^2, ...
       N11^1*N21^2*N31^2 + N12^1*N22^2*N32^2, ...
       N11^2*N21^2*N31^2 + N12^2*N22^2*N32^2, ...
       N11^3*N21^2*N31^2 + N12^3*N22^2*N32^2, ...
       N11^0*N21^3*N31^2 + N12^0*N22^3*N32^2, ...  
       N11^1*N21^3*N31^2 + N12^1*N22^3*N32^2, ...
       N11^2*N21^3*N31^2 + N12^2*N22^3*N32^2, ...
       N11^3*N21^3*N31^2 + N12^3*N22^3*N32^2, ...
       N11^0*N21^0*N31^3 + N12^0*N22^0*N32^3, ...
       N11^1*N21^0*N31^3 + N12^1*N22^0*N32^3, ...
       N11^2*N21^0*N31^3 + N12^2*N22^0*N32^3, ...
       N11^3*N21^0*N31^3 + N12^3*N22^0*N32^3, ...
       N11^0*N21^1*N31^3 + N12^0*N22^1*N32^3, ...
       N11^1*N21^1*N31^3 + N12^1*N22^1*N32^3, ...
       N11^2*N21^1*N31^3 + N12^2*N22^1*N32^3, ...
       N11^3*N21^1*N31^3 + N12^3*N22^1*N32^3, ...
       N11^0*N21^2*N31^3 + N12^0*N22^2*N32^3, ...
       N11^1*N21^2*N31^3 + N12^1*N22^2*N32^3, ...
       N11^2*N21^2*N31^3 + N12^2*N22^2*N32^3, ...
       N11^3*N21^2*N31^3 + N12^3*N22^2*N32^3, ...
       N11^0*N21^3*N31^3 + N12^0*N22^3*N32^3, ...  
       N11^1*N21^3*N31^3 + N12^1*N22^3*N32^3, ...
       N11^2*N21^3*N31^3 + N12^2*N22^3*N32^3, ...
       N11^3*N21^3*N31^3 + N12^3*N22^3*N32^3 ];
 
notions =  [     S100 S200 S300 S010 S110 S210 S310 S020 S120 S220 S320 S030 S130 S230 S330 ...
            S001 S101 S201 S301 S011 S111 S211 S311 S021 S121 S221 S321 S031 S131 S231 S331 ...
            S002 S102 S202 S302 S012 S112 S212 S312 S022 S122 S222 S322 S032 S132 S232 S332 ...
            S003 S103 S203 S303 S013 S113 S213 S313 S023 S123 S223 S323 S033 S133 S233 S333]; 
        
        
for i=1:length(notions)
    
    S = subs(S,sums(i),notions(i));
    
end

SS = sym('SS%d%d',[2,2]);
P = charpoly(SS);  

SS11 = S(1,1); SS12 = S(1,2); 
SS21 = S(2,1); SS22 = S(2,2); 

P = subs(P);

[cp1,tp1] = coeffs(P(1),[s1 s2 s3]); 
[cp2,tp2] = coeffs(P(2),[s1 s2 s3]); 
[cp3,tp3] = coeffs(P(3),[s1 s2 s3]);

cp1_000 = 1;
cp1_i00 = 1; 
cp1_0i0 = 1;
cp1_00i = 1;
cp1_ii0 = 1; 
cp1_i0i = 1; 
cp1_0ii = 1; 
cp1_iii = 1; 

cp2_000 = cp2(tp2==1);
cp2_i00 = cp2(tp2==s1); 
cp2_0i0 = cp2(tp2==s2);
cp2_00i = cp2(tp2==s3);
cp2_ii0 = cp2(tp2==s1*s2); 
cp2_i0i = cp2(tp2==s1*s3);
cp2_0ii = cp2(tp2==s2*s3);
cp2_iii = cp2(tp2==s1*s2*s3); 

cp3_000 = cp3(tp3==1);
cp3_i00 = cp3(tp3==s1^2); 
cp3_0i0 = cp3(tp3==s2^2);
cp3_00i = cp3(tp3==s3^2);
cp3_ii0 = cp3(tp3==s1^2*s2^2); 
cp3_i0i = cp3(tp3==s1^2*s3^2);
cp3_0ii = cp3(tp3==s2^2*s3^2);
cp3_iii = cp3(tp3==s1^2*s2^2*s3^2); 

S001 = simplify(Sigma_001); S010 = simplify(Sigma_010); S011 = simplify(Sigma_011); 
S100 = simplify(Sigma_100); S101 = simplify(Sigma_101); S110 = simplify(Sigma_110); 
S111 = simplify(Sigma_111); S200 = simplify(Sigma_200); S201 = simplify(Sigma_201); 
S210 = simplify(Sigma_210); S211 = simplify(Sigma_211); S300 = simplify(Sigma_300);
S301 = simplify(Sigma_301); S310 = simplify(Sigma_310); S311 = simplify(Sigma_311); 

F_cp2_000 = matlabFunction(subs(cp2_000));
F_cp2_i00 = matlabFunction(subs(cp2_i00));
F_cp2_0i0 = matlabFunction(subs(cp2_0i0));
F_cp2_ii0 = matlabFunction(subs(cp2_ii0));
F_cp2_00i = matlabFunction(subs(cp2_00i));
F_cp2_i0i = matlabFunction(subs(cp2_i0i));
F_cp2_0ii = matlabFunction(subs(cp2_0ii));
F_cp2_iii = matlabFunction(subs(cp2_iii));

F_cp3_000 = matlabFunction(subs(cp3_000));
F_cp3_i00 = matlabFunction(subs(cp3_i00));
F_cp3_0i0 = matlabFunction(subs(cp3_0i0));
F_cp3_ii0 = matlabFunction(subs(cp3_ii0));
F_cp3_00i = matlabFunction(subs(cp3_00i));
F_cp3_i0i = matlabFunction(subs(cp3_i0i));
F_cp3_0ii = matlabFunction(subs(cp3_0ii));
F_cp3_iii = matlabFunction(subs(cp3_iii));


%% Test 1

sample = 2^16; 

    a11 =  0.5; 
    a12 = -1.5; 
    a13 = -0.5; 
    a22 = 2.6; 
    a23 = -5; 
    a31 = -0.5; 
    a32 = -10; 
    a33 = 1; 
    r1 = 0.5; 
    r2 = -1.5; 
    r3 = -0.5; 
    b1 = 0.2; 
    b2 = -0.1; 
   
     a21 = linspace(-7,-1,sqrt(sample))';
     b3 = linspace(1.5,5,sqrt(sample))'; 
   
     [A,H] = meshgrid(a21,b3); 
     a21 = reshape(A,sample,1); 
     b3 = reshape(H,sample,1); 
     
    Cp2_000 = sign(F_cp2_000(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp2_i00 = sign(F_cp2_i00(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    Cp2_0i0 = sign(F_cp2_0i0(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    Cp2_ii0 = sign(F_cp2_ii0(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));  
    Cp2_00i = sign(F_cp2_00i(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp2_i0i = sign(F_cp2_i0i(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    Cp2_0ii = sign(F_cp2_0ii(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    Cp2_iii = sign(F_cp2_iii(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));  
    
    Cp3_000 = sign(F_cp3_000(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_i00 = sign(F_cp3_i00(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_0i0 = sign(F_cp3_0i0(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_ii0 = sign(F_cp3_ii0(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    Cp3_00i = sign(F_cp3_00i(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_i0i = sign(F_cp3_i0i(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_0ii = sign(F_cp3_0ii(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3)); 
    Cp3_iii = sign(F_cp3_iii(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3));
    
V000 = (1-sign(Cp2_000))/2 + (1-sign(Cp2_000).*sign(Cp3_000))/2;
Vi00 = (1-sign(Cp2_i00))/2 + (1-sign(Cp2_i00).*sign(Cp3_i00))/2; 
V0i0 = (1-sign(Cp2_0i0))/2 + (1-sign(Cp2_0i0).*sign(Cp3_0i0))/2; 
Vii0 = (1-sign(Cp2_ii0))/2 + (1-sign(Cp2_ii0).*sign(Cp3_ii0))/2; 
V00i = (1-sign(Cp2_00i))/2 + (1-sign(Cp2_00i).*sign(Cp3_00i))/2; 
Vi0i = (1-sign(Cp2_i0i))/2 + (1-sign(Cp2_i0i).*sign(Cp3_i0i))/2; 
V0ii = (1-sign(Cp2_0ii))/2 + (1-sign(Cp2_0ii).*sign(Cp3_0ii))/2; 
Viii = (1-sign(Cp2_iii))/2 + (1-sign(Cp2_iii).*sign(Cp3_iii))/2; 

NpRoots = (V000 - Vi00 - V0i0 - V00i + Vii0 + Vi0i + V0ii - Viii)/4; 

Signs = [Cp2_000, Cp2_i00, Cp2_0i0, Cp2_ii0, Cp2_00i, Cp2_i0i, Cp2_0ii, Cp2_iii, ... 
         Cp3_000, Cp3_i00, Cp3_0i0, Cp3_ii0, Cp3_00i, Cp3_i0i, Cp3_0ii, Cp3_iii, NpRoots];

unique_Signs = unique(Signs,'rows');
[rows,~] = find(unique_Signs(:,end)>0);
unique_nonzero_Signs = sortrows(unique_Signs(rows,1:end-1))

NPROOTS = transpose(reshape(NpRoots,sqrt(sample),sqrt(sample))); 
figure(1)
imagesc(b3,a21,NPROOTS);
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('b_{3}', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Number of feasible roots', 'FontSize', 30); 
annotation('textbox', [0.6, 0.25, 0.1, 0.1], 'String', ["Orange region: #of feasible roots = 1", "White region: #of feasible roots = 0"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

figure(2)
Xf=F_cp2_i00(a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,r1,r2,r3);
Xf=(Xf>0 ==1) & (Xf<0 == 0);
XF = transpose(reshape(Xf,sqrt(sample),sqrt(sample))); 
imagesc(b3,a21,XF); colorbar off
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('b_{3}', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Sign of v_1(\infty,0,0)', 'FontSize', 30); 
annotation('textbox', [0.6, 0.25, 0.1, 0.1], 'String', ["Orange region: positive sign", "White region: negative sign"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')
unique(XF-NPROOTS)


%%
clc

    a11 =  0.5; 
    a12 = -1.5; 
    a13 = -0.5; 
    a22 = 2.6; 
    a23 = -5; 
    a31 = -0.5; 
    a32 = -10; 
    a33 = 1; 
    r1 = 0.5; 
    r2 = -1.5; 
    r3 = -0.5; 
    b1 = 0.2; 
    b2 = -0.1;   
    
    a21 = -1.2;
    b3 = 4; 

set_phcpath('/Users/mohammadaladwani/Documents/MATLAB/PHClab1.0.4/phc'); 

MM = [r1,  0, 0, 0; ...
      a11, 1, 0, 0; ...
      a12, 0, 1, 0; ...
      a13, 0, 0, 1; ...
      b1,  1, 1, 0; ...
      0,   0, 0, 0; ...
      r2,  0, 0, 0; ...
      a21, 1, 0, 0; ...
      a22, 0, 1, 0; ...
      a23, 0, 0, 1; ...
      b2 , 1, 1, 0; ...
      0,   0, 0, 0; ...
      r3,  0, 0, 0; ...
      a31, 1, 0, 0; ...
      a32, 0, 1, 0; ...
      a33, 0, 0, 1; ...
      b3 , 1, 1, 0; ...
      0,   0, 0, 0];    
  
make_system(MM); % shows symbolic format of the system
s = solve_system(MM); % call the blackbox solver
ns = size(s,2); % check the number of solutions
 
X = [s.x1; s.x2; s.x3]; 
count_feasible = 0;
size(X,2);

for i = 1:size(X,2)
    
    if (real(X(1,i)) > 10^-6 && real(X(2,i)) > 10^-6 && real(X(3,i)) > 10^-6 && ... 
    abs(imag(X(1,i))) < 10^-6 && abs(imag(X(2,i))) < 10^-6 && abs(imag(X(3,i))) < 10^-6)

         count_feasible  = count_feasible  + 1;
    
    end

    
end

count_feasible

check_sum = 0;
k1 = 1;
k2 = 1; 
k3 = 0;

for i = 1:size(X,2)
    
    Y1 = X(1,i);
    Y2 = X(2,i);
    Y3 = X(3,i); 

         check_sum  = check_sum  + Y1^k1*Y2^k2*Y3^k3;
    
end

%check_S222_1 = vpa(subs((sprintf('Sigma_%d%d%d',k1,k2,k3))),6)
check_S222_2 = vpa(check_sum,6)


