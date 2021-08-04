clear all; close all; clc; 

tic
digits(100); 

    a11 =  2; 
    a12 = -1.5; 
    a13 = -1.5; 
    a22 = 2; 
    a23 = -1.5; 
    a31 = -1.5; 
    a32 = -1; 
    a33 = 1; 
    r1 = 1.5; 
    r2 = -1.5; 
    r3 = -1.5; 
    b1 = 1; 
    b2 = -1; 
    
% part 1

syms N1 N2 N3 ...
     F1 F2 F3 ... 
     a21 b3 ...
 
f1 = @(N1,N2,N3) r1+a11*N1+a12*N2+a13*N3+b1*N2*N3;
f2 = @(N1,N2,N3) r2+a21*N1+a22*N2+a23*N3+b2*N1*N3;
f3 = @(N1,N2,N3) r3+a31*N1+a32*N2+a33*N3+b3*N1*N2;

r = 6; 

HrowL1 = [f1(N1,N2,N3), f2(N1,N2,N3), f2(N1,N2,N3)*N2, f2(N1,N2,N3)*N3,...
          f3(N1,N2,N3), f3(N1,N2,N3)*N3]; 
      
HcolL1 = [1, N2, N2^2, N2*N3, N3, N3^2];  

for i = 1:r
    
    [c,T] = coeffs(HrowL1(i),[N2,N3]);
    
    for j = 1:r
        
        if isempty(c(T==HcolL1(j))) 
            
            M_N1(i,j) = 0;
            
        else
            
            M_N1(i,j) = c(T==HcolL1(j)); 
            
        end
        
    end
    
end

Res_N2N3 = collect(det(M_N1),N1);

Mp_N1 = M_N1;
Mp_N1(:,1) = [F1, F2, F2*N2, F2*N3, F3, F3*N3]; 

T1 = collect(det(Mp_N1),[F1 F2 F3]);
T11 = subs(T1, [F1 F2 F3], [1 0 0]); 
T12 = subs(T1, [F1 F2 F3], [0 1 0]); 
T13 = subs(T1, [F1 F2 F3], [0 0 1]); 

t11 = matlabFunction(T11); 
fT11 =@(N1)       t11(N1,b3);
t12 = matlabFunction(T12); 
fT12 =@(N1,N3)    t12(N1,N3,a21,b3);
t13 = matlabFunction(T13); 
fT13 =@(N1,N3)    t13(N1,N3,a21,b3);

%
HrowL2 = [f1(N1,N2,N3), f1(N1,N2,N3)*N1, f1(N1,N2,N3)*N3, f2(N1,N2,N3),...
          f3(N1,N2,N3)*N1, f3(N1,N2,N3)*N3]; 
      
HcolL2 = [1, N1, N3, N1^2, N1*N3, N3^2];  

for i = 1:r
    
    [c,T] = coeffs(HrowL2(i),[N1,N3]);
    
    for j = 1:r
        
        if isempty(c(T==HcolL2(j))) 
            
            M_N2(i,j) = 0;
            
        else
            
            M_N2(i,j) = c(T==HcolL2(j)); 
            
        end
        
    end
    
end

Res_N1N3 = collect(det(M_N2),N2);

Mp_N2 = M_N2;
Mp_N2(:,1) = [F1, F1*N1, F1*N3, F2, F3*N1, F3*N3];

T2 = collect(det(Mp_N2),[F1 F2 F3]);
T21 = subs(T2, [F1 F2 F3], [1 0 0]); 
T22 = subs(T2, [F1 F2 F3], [0 1 0]); 
T23 = subs(T2, [F1 F2 F3], [0 0 1]); 

t21 = matlabFunction(T21); 
fT21 =@(N1,N2,N3) t21(N1,N2,N3,a21,b3);
t22 = matlabFunction(T22); 
fT22 =@(N2)       t22(N2,b3);
t23 = matlabFunction(T23); 
fT23 =@(N1,N2,N3) t23(N1,N2,N3,a21,b3);

%
HrowL3 = [f1(N1,N2,N3), f1(N1,N2,N3)*N1, f1(N1,N2,N3)*N2, f2(N1,N2,N3)*N1,...
          f2(N1,N2,N3)*N2, f3(N1,N2,N3)]; 
      
HcolL3 = [1, N1, N2, N1^2, N1*N2, N2^2];  

for i = 1:r
    
    [c,T] = coeffs(HrowL3(i),[N1,N2]);
    
    for j = 1:r
        
        if isempty(c(T==HcolL3(j))) 
            
            M_N3(i,j) = 0;
            
        else
            
            M_N3(i,j) = c(T==HcolL3(j)); 
            
        end
        
    end
    
end

Res_N1N2 = collect(det(M_N3),N3);

Mp_N3 = M_N3;
Mp_N3(:,1) = [F1, F1*N1, F1*N2, F2*N1, F2*N2, F3];

T3 = collect(det(Mp_N3),[F1 F2 F3]);
T31 = subs(T3, [F1 F2 F3], [1 0 0]); 
T32 = subs(T3, [F1 F2 F3], [0 1 0]); 
T33 = subs(T3, [F1 F2 F3], [0 0 1]); 

t31 = matlabFunction(T31); 
fT31 =@(N1,N2,N3) t31(N1,N2,N3,a21,b3);
t32 = matlabFunction(T32); 
fT32 =@(N1,N2,N3) t32(N1,N2,N3,a21,b3);
t33 = matlabFunction(T33); 
fT33 =@(N3)       t33(N3,a21);


T =@(N1,N2,N3) det([fT11(N1)       fT12(N1,N3)    fT13(N1,N3); ...
                    fT21(N1,N2,N3) fT22(N2)       fT23(N1,N2,N3); ...
                    fT31(N1,N2,N3) fT32(N1,N2,N3) fT33(N3)]); 
                
Jacob = det(jacobian([f1(N1,N2,N3),f2(N1,N2,N3),f3(N1,N2,N3)],[N1,N2,N3]));
fJ = matlabFunction(Jacob);
J =@(N1,N2,N3) fJ(N1,N2,N3,a21,b3);

R1 = matlabFunction(Res_N2N3);
Res1 =@(N1) R1(N1,a21,b3);
R2 = matlabFunction(Res_N1N3);
Res2 =@(N2) R2(N2,a21,b3);
R3 = matlabFunction(Res_N1N2);
Res3 =@(N3) R3(N3,a21,b3);

syms x y z

H1 = 1/Res1(1/x);
L1 = simplify(taylor(H1,'Order', 15));
LL1 = matlabFunction(L1); 
Resultant1 =@(x) LL1(a21,b3,x);

H2 = 1/Res2(1/y);
L2 = simplify(taylor(H2,'Order', 15));
LL2 = matlabFunction(L2); 
Resultant2 =@(y) LL2(a21,b3,y);

H3 = 1/Res3(1/z);
L3 = simplify(taylor(H3,'Order', 15));
LL3 = matlabFunction(L3); 
Resultant3 =@(z) LL3(a21,b3,z);
    
[cTJ, tTJ] = coeffs(T(x,y,z)*J(x,y,z),[x y z]);
[cx, tx] = coeffs(Resultant1(x),x);
[cy, ty] = coeffs(Resultant2(y),y);
[cz, tz] = coeffs(Resultant3(z),z);

Sigma_000 = 0; Sigma_100 = 0; Sigma_200 = 0; Sigma_300 = 0;
Sigma_010 = 0; Sigma_110 = 0; Sigma_210 = 0; Sigma_310 = 0;
Sigma_020 = 0; Sigma_120 = 0; Sigma_220 = 0; Sigma_320 = 0;
Sigma_030 = 0; Sigma_130 = 0; Sigma_230 = 0; Sigma_330 = 0;

Sigma_001 = 0; Sigma_101 = 0; Sigma_201 = 0; Sigma_301 = 0;
Sigma_011 = 0; Sigma_111 = 0; Sigma_211 = 0; Sigma_311 = 0;
Sigma_021 = 0; Sigma_121 = 0; Sigma_221 = 0; Sigma_321 = 0;
Sigma_031 = 0; Sigma_131 = 0; Sigma_231 = 0; Sigma_331 = 0;

Sigma_002 = 0; Sigma_102 = 0; Sigma_202 = 0; Sigma_302 = 0;
Sigma_012 = 0; Sigma_112 = 0; Sigma_212 = 0; Sigma_312 = 0;
Sigma_022 = 0; Sigma_122 = 0; Sigma_222 = 0; Sigma_322 = 0;
Sigma_032 = 0; Sigma_132 = 0; Sigma_232 = 0; Sigma_332 = 0;
                                            
Sigma_003 = 0; Sigma_103 = 0; Sigma_203 = 0; Sigma_303 = 0;
Sigma_013 = 0; Sigma_113 = 0; Sigma_213 = 0; Sigma_313 = 0;
Sigma_023 = 0; Sigma_123 = 0; Sigma_223 = 0; Sigma_323 = 0;
Sigma_033 = 0; Sigma_133 = 0; Sigma_233 = 0; Sigma_333 = 0;              
                                             
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
    
        
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_020 = Sigma_020 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_120 = Sigma_120 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_220 = Sigma_220 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_320 = Sigma_320 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+0)); 
    end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_030 = Sigma_030 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_130 = Sigma_130 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_230 = Sigma_230 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+0)))==1
        Sigma_330 = Sigma_330 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+0)); 
    end
    
    %
    
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
    
        
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_021 = Sigma_021 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_121 = Sigma_121 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_221 = Sigma_221 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_321 = Sigma_321 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+1)); 
    end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_031 = Sigma_031 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_131 = Sigma_131 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_231 = Sigma_231 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+1)))==1
        Sigma_331 = Sigma_331 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+1)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_002 = Sigma_002 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_102 = Sigma_102 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_202 = Sigma_202 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+2)); 
    end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_302 = Sigma_302 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+2)); 
     end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_012 = Sigma_012 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_112 = Sigma_112 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_212 = Sigma_212 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+2)); 
    end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_312 = Sigma_312 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+2)); 
     end
    
        
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_022 = Sigma_022 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_122 = Sigma_122 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+2)))==1
        Sigma_222 = Sigma_222 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+2)); 
    end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_322 = Sigma_322 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+2)); 
     end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_032 = Sigma_032 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+2)); 
    end
    
     if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_132 = Sigma_132 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+2)); 
     end
    
     if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_232 = Sigma_232 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+2)); 
     end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+2)))==1
         Sigma_332 = Sigma_332 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+2)); 
     end

    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+3)))==1
        Sigma_003 = Sigma_003 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+3)))==1
        Sigma_103 = Sigma_103 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+3)); 
    end
    
     if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_203 = Sigma_203 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+3)); 
     end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+0)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_303 = Sigma_303 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+0))*cz(tz==z^(powz+1+3)); 
     end
    
    
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+3)))==1
        Sigma_013 = Sigma_013 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+3)))==1
        Sigma_113 = Sigma_113 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+3)); 
    end
    
     if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_213 = Sigma_213 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+3)); 
     end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+1)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_313 = Sigma_313 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+1))*cz(tz==z^(powz+1+3)); 
     end
    
        
    
     if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_023 = Sigma_023 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+3)); 
     end
   
     if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_123 = Sigma_123 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+3)); 
     end
    
      if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_223 = Sigma_223 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+3)); 
      end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+2)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_323 = Sigma_323 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+2))*cz(tz==z^(powz+1+3)); 
     end
    
    

     if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_033 = Sigma_033 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+3)); 
     end
    
     if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_133 = Sigma_133 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+3)); 
     end
    
     if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_233 = Sigma_233 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+3)); 
     end
    
     if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+3)))*max(ismember(tz,z^(powz+1+3)))==1
         Sigma_333 = Sigma_333 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+3))*cz(tz==z^(powz+1+3)); 
     end

    
end

% part 2

syms N11 N12 N13 N14 N15 N21 N22 N23 N24 N25 N31 N32 N33 N34 N35 s1 s2 s3

%m=@(N1,N2,N3) [1; N1; N2; N3; N1*N2];
m=@(N1,N2,N3) [1; N1; N1*N2; N1*N3; N1*N2*N3];

W=[m(N11,N21,N31),m(N12,N22,N32),m(N13,N23,N33),m(N14,N24,N34),m(N15,N25,N35)]; 
D = diag([(N11-s1)*(N21-s2)*(N31-s3),(N12-s1)*(N22-s2)*(N32-s3),(N13-s1)*(N23-s2)*(N33-s3),(N14-s1)*(N24-s2)*(N34-s3),(N15-s1)*(N25-s2)*(N35-s3)]);
S=collect(W*D*transpose(W),[s1 s2 s3]); 

syms  S000 S100 S200 S300 S010 S110 S210 S310 S020 S120 S220 S320 S030 S130 S230 S330 ...
      S001 S101 S201 S301 S011 S111 S211 S311 S021 S121 S221 S321 S031 S131 S231 S331 ...
      S002 S102 S202 S302 S012 S112 S212 S312 S022 S122 S222 S322 S032 S132 S232 S332 ...
      S003 S103 S203 S303 S013 S113 S213 S313 S023 S123 S223 S323 S033 S133 S233 S333
  
sums =[N11^1*N21^0*N31^0 + N12^1*N22^0*N32^0 + N13^1*N23^0*N33^0 + N14^1*N24^0*N34^0 + N15^1*N25^0*N35^0, ...
       N11^2*N21^0*N31^0 + N12^2*N22^0*N32^0 + N13^2*N23^0*N33^0 + N14^2*N24^0*N34^0 + N15^2*N25^0*N35^0, ...
       N11^3*N21^0*N31^0 + N12^3*N22^0*N32^0 + N13^3*N23^0*N33^0 + N14^3*N24^0*N34^0 + N15^3*N25^0*N35^0, ...
       N11^0*N21^1*N31^0 + N12^0*N22^1*N32^0 + N13^0*N23^1*N33^0 + N14^0*N24^1*N34^0 + N15^0*N25^1*N35^0, ...
       N11^1*N21^1*N31^0 + N12^1*N22^1*N32^0 + N13^1*N23^1*N33^0 + N14^1*N24^1*N34^0 + N15^1*N25^1*N35^0, ...
       N11^2*N21^1*N31^0 + N12^2*N22^1*N32^0 + N13^2*N23^1*N33^0 + N14^2*N24^1*N34^0 + N15^2*N25^1*N35^0, ...
       N11^3*N21^1*N31^0 + N12^3*N22^1*N32^0 + N13^3*N23^1*N33^0 + N14^3*N24^1*N34^0 + N15^3*N25^1*N35^0, ...
       N11^0*N21^2*N31^0 + N12^0*N22^2*N32^0 + N13^0*N23^2*N33^0 + N14^0*N24^2*N34^0 + N15^0*N25^2*N35^0, ...
       N11^1*N21^2*N31^0 + N12^1*N22^2*N32^0 + N13^1*N23^2*N33^0 + N14^1*N24^2*N34^0 + N15^1*N25^2*N35^0, ...
       N11^2*N21^2*N31^0 + N12^2*N22^2*N32^0 + N13^2*N23^2*N33^0 + N14^2*N24^2*N34^0 + N15^2*N25^2*N35^0, ...
       N11^3*N21^2*N31^0 + N12^3*N22^2*N32^0 + N13^3*N23^2*N33^0 + N14^3*N24^2*N34^0 + N15^3*N25^2*N35^0, ...
       N11^0*N21^3*N31^0 + N12^0*N22^3*N32^0 + N13^0*N23^3*N33^0 + N14^0*N24^3*N34^0 + N15^0*N25^3*N35^0, ...  
       N11^1*N21^3*N31^0 + N12^1*N22^3*N32^0 + N13^1*N23^3*N33^0 + N14^1*N24^3*N34^0 + N15^1*N25^3*N35^0, ...
       N11^2*N21^3*N31^0 + N12^2*N22^3*N32^0 + N13^2*N23^3*N33^0 + N14^2*N24^3*N34^0 + N15^2*N25^3*N35^0, ...
       N11^3*N21^3*N31^0 + N12^3*N22^3*N32^0 + N13^3*N23^3*N33^0 + N14^3*N24^3*N34^0 + N15^3*N25^3*N35^0, ...   
       N11^0*N21^0*N31^1 + N12^0*N22^0*N32^1 + N13^0*N23^0*N33^1 + N14^0*N24^0*N34^1 + N15^0*N25^0*N35^1, ...
       N11^1*N21^0*N31^1 + N12^1*N22^0*N32^1 + N13^1*N23^0*N33^1 + N14^1*N24^0*N34^1 + N15^1*N25^0*N35^1, ...
       N11^2*N21^0*N31^1 + N12^2*N22^0*N32^1 + N13^2*N23^0*N33^1 + N14^2*N24^0*N34^1 + N15^2*N25^0*N35^1, ...
       N11^3*N21^0*N31^1 + N12^3*N22^0*N32^1 + N13^3*N23^0*N33^1 + N14^3*N24^0*N34^1 + N15^3*N25^0*N35^1, ...
       N11^0*N21^1*N31^1 + N12^0*N22^1*N32^1 + N13^0*N23^1*N33^1 + N14^0*N24^1*N34^1 + N15^0*N25^1*N35^1, ...
       N11^1*N21^1*N31^1 + N12^1*N22^1*N32^1 + N13^1*N23^1*N33^1 + N14^1*N24^1*N34^1 + N15^1*N25^1*N35^1, ...
       N11^2*N21^1*N31^1 + N12^2*N22^1*N32^1 + N13^2*N23^1*N33^1 + N14^2*N24^1*N34^1 + N15^2*N25^1*N35^1, ...
       N11^3*N21^1*N31^1 + N12^3*N22^1*N32^1 + N13^3*N23^1*N33^1 + N14^3*N24^1*N34^1 + N15^3*N25^1*N35^1, ...
       N11^0*N21^2*N31^1 + N12^0*N22^2*N32^1 + N13^0*N23^2*N33^1 + N14^0*N24^2*N34^1 + N15^0*N25^2*N35^1, ...
       N11^1*N21^2*N31^1 + N12^1*N22^2*N32^1 + N13^1*N23^2*N33^1 + N14^1*N24^2*N34^1 + N15^1*N25^2*N35^1, ...
       N11^2*N21^2*N31^1 + N12^2*N22^2*N32^1 + N13^2*N23^2*N33^1 + N14^2*N24^2*N34^1 + N15^2*N25^2*N35^1, ...
       N11^3*N21^2*N31^1 + N12^3*N22^2*N32^1 + N13^3*N23^2*N33^1 + N14^3*N24^2*N34^1 + N15^3*N25^2*N35^1, ...
       N11^0*N21^3*N31^1 + N12^0*N22^3*N32^1 + N13^0*N23^3*N33^1 + N14^0*N24^3*N34^1 + N15^0*N25^3*N35^1, ...  
       N11^1*N21^3*N31^1 + N12^1*N22^3*N32^1 + N13^1*N23^3*N33^1 + N14^1*N24^3*N34^1 + N15^1*N25^3*N35^1, ...
       N11^2*N21^3*N31^1 + N12^2*N22^3*N32^1 + N13^2*N23^3*N33^1 + N14^2*N24^3*N34^1 + N15^2*N25^3*N35^1, ...
       N11^3*N21^3*N31^1 + N12^3*N22^3*N32^1 + N13^3*N23^3*N33^1 + N14^3*N24^3*N34^1 + N15^3*N25^3*N35^1, ...
       N11^0*N21^0*N31^2 + N12^0*N22^0*N32^2 + N13^0*N23^0*N33^2 + N14^0*N24^0*N34^2 + N15^0*N25^0*N35^2, ...
       N11^1*N21^0*N31^2 + N12^1*N22^0*N32^2 + N13^1*N23^0*N33^2 + N14^1*N24^0*N34^2 + N15^1*N25^0*N35^2, ...
       N11^2*N21^0*N31^2 + N12^2*N22^0*N32^2 + N13^2*N23^0*N33^2 + N14^2*N24^0*N34^2 + N15^2*N25^0*N35^2, ...
       N11^3*N21^0*N31^2 + N12^3*N22^0*N32^2 + N13^3*N23^0*N33^2 + N14^3*N24^0*N34^2 + N15^3*N25^0*N35^2, ...
       N11^0*N21^1*N31^2 + N12^0*N22^1*N32^2 + N13^0*N23^1*N33^2 + N14^0*N24^1*N34^2 + N15^0*N25^1*N35^2, ...
       N11^1*N21^1*N31^2 + N12^1*N22^1*N32^2 + N13^1*N23^1*N33^2 + N14^1*N24^1*N34^2 + N15^1*N25^1*N35^2, ...
       N11^2*N21^1*N31^2 + N12^2*N22^1*N32^2 + N13^2*N23^1*N33^2 + N14^2*N24^1*N34^2 + N15^2*N25^1*N35^2, ...
       N11^3*N21^1*N31^2 + N12^3*N22^1*N32^2 + N13^3*N23^1*N33^2 + N14^3*N24^1*N34^2 + N15^3*N25^1*N35^2, ...
       N11^0*N21^2*N31^2 + N12^0*N22^2*N32^2 + N13^0*N23^2*N33^2 + N14^0*N24^2*N34^2 + N15^0*N25^2*N35^2, ...
       N11^1*N21^2*N31^2 + N12^1*N22^2*N32^2 + N13^1*N23^2*N33^2 + N14^1*N24^2*N34^2 + N15^1*N25^2*N35^2, ...
       N11^2*N21^2*N31^2 + N12^2*N22^2*N32^2 + N13^2*N23^2*N33^2 + N14^2*N24^2*N34^2 + N15^2*N25^2*N35^2, ...
       N11^3*N21^2*N31^2 + N12^3*N22^2*N32^2 + N13^3*N23^2*N33^2 + N14^3*N24^2*N34^2 + N15^3*N25^2*N35^2, ...
       N11^0*N21^3*N31^2 + N12^0*N22^3*N32^2 + N13^0*N23^3*N33^2 + N14^0*N24^3*N34^2 + N15^0*N25^3*N35^2, ...  
       N11^1*N21^3*N31^2 + N12^1*N22^3*N32^2 + N13^1*N23^3*N33^2 + N14^1*N24^3*N34^2 + N15^1*N25^3*N35^2, ...
       N11^2*N21^3*N31^2 + N12^2*N22^3*N32^2 + N13^2*N23^3*N33^2 + N14^2*N24^3*N34^2 + N15^2*N25^3*N35^2, ...
       N11^3*N21^3*N31^2 + N12^3*N22^3*N32^2 + N13^3*N23^3*N33^2 + N14^3*N24^3*N34^2 + N15^3*N25^3*N35^2, ...
       N11^0*N21^0*N31^3 + N12^0*N22^0*N32^3 + N13^0*N23^0*N33^3 + N14^0*N24^0*N34^3 + N15^0*N25^0*N35^3, ...
       N11^1*N21^0*N31^3 + N12^1*N22^0*N32^3 + N13^1*N23^0*N33^3 + N14^1*N24^0*N34^3 + N15^1*N25^0*N35^3, ...
       N11^2*N21^0*N31^3 + N12^2*N22^0*N32^3 + N13^2*N23^0*N33^3 + N14^2*N24^0*N34^3 + N15^2*N25^0*N35^3, ...
       N11^3*N21^0*N31^3 + N12^3*N22^0*N32^3 + N13^3*N23^0*N33^3 + N14^3*N24^0*N34^3 + N15^3*N25^0*N35^3, ...
       N11^0*N21^1*N31^3 + N12^0*N22^1*N32^3 + N13^0*N23^1*N33^3 + N14^0*N24^1*N34^3 + N15^0*N25^1*N35^3, ...
       N11^1*N21^1*N31^3 + N12^1*N22^1*N32^3 + N13^1*N23^1*N33^3 + N14^1*N24^1*N34^3 + N15^1*N25^1*N35^3, ...
       N11^2*N21^1*N31^3 + N12^2*N22^1*N32^3 + N13^2*N23^1*N33^3 + N14^2*N24^1*N34^3 + N15^2*N25^1*N35^3, ...
       N11^3*N21^1*N31^3 + N12^3*N22^1*N32^3 + N13^3*N23^1*N33^3 + N14^3*N24^1*N34^3 + N15^3*N25^1*N35^3, ...
       N11^0*N21^2*N31^3 + N12^0*N22^2*N32^3 + N13^0*N23^2*N33^3 + N14^0*N24^2*N34^3 + N15^0*N25^2*N35^3, ...
       N11^1*N21^2*N31^3 + N12^1*N22^2*N32^3 + N13^1*N23^2*N33^3 + N14^1*N24^2*N34^3 + N15^1*N25^2*N35^3, ...
       N11^2*N21^2*N31^3 + N12^2*N22^2*N32^3 + N13^2*N23^2*N33^3 + N14^2*N24^2*N34^3 + N15^2*N25^2*N35^3, ...
       N11^3*N21^2*N31^3 + N12^3*N22^2*N32^3 + N13^3*N23^2*N33^3 + N14^3*N24^2*N34^3 + N15^3*N25^2*N35^3, ...
       N11^0*N21^3*N31^3 + N12^0*N22^3*N32^3 + N13^0*N23^3*N33^3 + N14^0*N24^3*N34^3 + N15^0*N25^3*N35^3, ...  
       N11^1*N21^3*N31^3 + N12^1*N22^3*N32^3 + N13^1*N23^3*N33^3 + N14^1*N24^3*N34^3 + N15^1*N25^3*N35^3, ...
       N11^2*N21^3*N31^3 + N12^2*N22^3*N32^3 + N13^2*N23^3*N33^3 + N14^2*N24^3*N34^3 + N15^2*N25^3*N35^3, ...
       N11^3*N21^3*N31^3 + N12^3*N22^3*N32^3 + N13^3*N23^3*N33^3 + N14^3*N24^3*N34^3 + N15^3*N25^3*N35^3];
 
notions =  [     S100 S200 S300 S010 S110 S210 S310 S020 S120 S220 S320 S030 S130 S230 S330 ...
            S001 S101 S201 S301 S011 S111 S211 S311 S021 S121 S221 S321 S031 S131 S231 S331 ...
            S002 S102 S202 S302 S012 S112 S212 S312 S022 S122 S222 S322 S032 S132 S232 S332 ...
            S003 S103 S203 S303 S013 S113 S213 S313 S023 S123 S223 S323 S033 S133 S233 S333]; 
        
        
for i=1:length(notions)
    
    S = subs(S,sums(i),notions(i));
    
end

SS = sym('SS%d%d',[5,5]);
P = charpoly(SS);  

SS11 = S(1,1); SS12 = S(1,2); SS13 = S(1,3); SS14 = S(1,4); SS15 = S(1,5); 
SS21 = S(2,1); SS22 = S(2,2); SS23 = S(2,3); SS24 = S(2,4); SS25 = S(2,5); 
SS31 = S(3,1); SS32 = S(3,2); SS33 = S(3,3); SS34 = S(3,4); SS35 = S(3,5); 
SS41 = S(4,1); SS42 = S(4,2); SS43 = S(4,3); SS44 = S(4,4); SS45 = S(4,5);
SS51 = S(5,1); SS52 = S(5,2); SS53 = S(5,3); SS54 = S(5,4); SS55 = S(5,5); 

P = subs(P);

[cp1,tp1] = coeffs(P(1),[s1 s2 s3]); 
[cp2,tp2] = coeffs(P(2),[s1 s2 s3]); 
[cp3,tp3] = coeffs(P(3),[s1 s2 s3]);
[cp4,tp4] = coeffs(P(4),[s1 s2 s3]); 
[cp5,tp5] = coeffs(P(5),[s1 s2 s3]);
[cp6,tp6] = coeffs(P(6),[s1 s2 s3]);

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

cp4_000 = cp4(tp4==1);
cp4_i00 = cp4(tp4==s1^3); 
cp4_0i0 = cp4(tp4==s2^3);
cp4_00i = cp4(tp4==s3^3);
cp4_ii0 = cp4(tp4==s1^3*s2^3); 
cp4_i0i = cp4(tp4==s1^3*s3^3);
cp4_0ii = cp4(tp4==s2^3*s3^3);
cp4_iii = cp4(tp4==s1^3*s2^3*s3^3); 

cp5_000 = cp5(tp5==1);
cp5_i00 = cp5(tp5==s1^4); 
cp5_0i0 = cp5(tp5==s2^4);
cp5_00i = cp5(tp5==s3^4);
cp5_ii0 = cp5(tp5==s1^4*s2^4); 
cp5_i0i = cp5(tp5==s1^4*s3^4);
cp5_0ii = cp5(tp5==s2^4*s3^4);
cp5_iii = cp5(tp5==s1^4*s2^4*s3^4); 

cp6_000 = cp6(tp6==1);
cp6_i00 = cp6(tp6==s1^5); 
cp6_0i0 = cp6(tp6==s2^5);
cp6_00i = cp6(tp6==s3^5);
cp6_ii0 = cp6(tp6==s1^5*s2^5); 
cp6_i0i = cp6(tp6==s1^5*s3^5);
cp6_0ii = cp6(tp6==s2^5*s3^5);
cp6_iii = cp6(tp6==s1^5*s2^5*s3^5); 

S000 = simplify(Sigma_000) ; S100 = simplify(Sigma_100) ; S200 = simplify(Sigma_200) ; S300 = simplify(Sigma_300) ;
S010 = simplify(Sigma_010) ; S110 = simplify(Sigma_110) ; S210 = simplify(Sigma_210) ; S310 = simplify(Sigma_310) ;
S020 = simplify(Sigma_020) ; S120 = simplify(Sigma_120) ; S220 = simplify(Sigma_220) ; S320 = simplify(Sigma_320) ;
S030 = simplify(Sigma_030) ; S130 = simplify(Sigma_130) ; S230 = simplify(Sigma_230) ; S330 = simplify(Sigma_330) ;

S001 = simplify(Sigma_001) ; S101 = simplify(Sigma_101) ; S201 = simplify(Sigma_201) ; S301 = simplify(Sigma_301) ;
S011 = simplify(Sigma_011) ; S111 = simplify(Sigma_111) ; S211 = simplify(Sigma_211) ; S311 = simplify(Sigma_311) ;
S021 = simplify(Sigma_021) ; S121 = simplify(Sigma_121) ; S221 = simplify(Sigma_221) ; S321 = simplify(Sigma_321) ;
S031 = simplify(Sigma_031) ; S131 = simplify(Sigma_131) ; S231 = simplify(Sigma_231) ; S331 = simplify(Sigma_331) ;

S002 = simplify(Sigma_002) ; S102 = simplify(Sigma_102) ; S202 = simplify(Sigma_202) ; S302 = simplify(Sigma_302) ;
S012 = simplify(Sigma_012) ; S112 = simplify(Sigma_112) ; S212 = simplify(Sigma_212) ; S312 = simplify(Sigma_312) ;
S022 = simplify(Sigma_022) ; S122 = simplify(Sigma_122) ; S222 = simplify(Sigma_222) ; S322 = simplify(Sigma_322) ;
S032 = simplify(Sigma_032) ; S132 = simplify(Sigma_132) ; S232 = simplify(Sigma_232) ; S332 = simplify(Sigma_332) ;
                                                                          
S003 = simplify(Sigma_003) ; S103 = simplify(Sigma_103) ; S203 = simplify(Sigma_203) ; S303 = simplify(Sigma_303) ;
S013 = simplify(Sigma_013) ; S113 = simplify(Sigma_113) ; S213 = simplify(Sigma_213) ; S313 = simplify(Sigma_313) ;
S023 = simplify(Sigma_023) ; S123 = simplify(Sigma_123) ; S223 = simplify(Sigma_223) ; S323 = simplify(Sigma_323) ;
S033 = simplify(Sigma_033) ; S133 = simplify(Sigma_133) ; S233 = simplify(Sigma_233) ; S333 = simplify(Sigma_333) ;

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

F_cp4_000 = matlabFunction(subs(cp4_000));
F_cp4_i00 = matlabFunction(subs(cp4_i00));
F_cp4_0i0 = matlabFunction(subs(cp4_0i0));
F_cp4_ii0 = matlabFunction(subs(cp4_ii0));
F_cp4_00i = matlabFunction(subs(cp4_00i));
F_cp4_i0i = matlabFunction(subs(cp4_i0i));
F_cp4_0ii = matlabFunction(subs(cp4_0ii));
F_cp4_iii = matlabFunction(subs(cp4_iii));

F_cp5_000 = matlabFunction(subs(cp5_000));
F_cp5_i00 = matlabFunction(subs(cp5_i00));
F_cp5_0i0 = matlabFunction(subs(cp5_0i0));
F_cp5_ii0 = matlabFunction(subs(cp5_ii0));
F_cp5_00i = matlabFunction(subs(cp5_00i));
F_cp5_i0i = matlabFunction(subs(cp5_i0i));
F_cp5_0ii = matlabFunction(subs(cp5_0ii));
F_cp5_iii = matlabFunction(subs(cp5_iii));

F_cp6_000 = matlabFunction(subs(cp6_000));
F_cp6_i00 = matlabFunction(subs(cp6_i00));
F_cp6_0i0 = matlabFunction(subs(cp6_0i0));
F_cp6_ii0 = matlabFunction(subs(cp6_ii0));
F_cp6_00i = matlabFunction(subs(cp6_00i));
F_cp6_i0i = matlabFunction(subs(cp6_i0i));
F_cp6_0ii = matlabFunction(subs(cp6_0ii));
F_cp6_iii = matlabFunction(subs(cp6_iii));

%
% Test

sample = 2^16; 
   
     a21 = linspace(1,6,sqrt(sample))';
     b3 = linspace(2,5,sqrt(sample))'; 
   
     [A,H] = meshgrid(a21,b3); 
     a21 = reshape(A,sample,1); 
     b3 = reshape(H,sample,1); 
     
    Cp2_000 = sign(F_cp2_000(a21,b3)); 
    Cp2_i00 = sign(F_cp2_i00(a21,b3));
    Cp2_0i0 = sign(F_cp2_0i0(a21,b3));
    Cp2_ii0 = sign(F_cp2_ii0(a21,b3));  
    Cp2_00i = sign(F_cp2_00i(a21,b3)); 
    Cp2_i0i = sign(F_cp2_i0i(a21,b3));
    Cp2_0ii = sign(F_cp2_0ii(a21,b3));
    Cp2_iii = sign(F_cp2_iii(a21,b3));  
    
    Cp3_000 = sign(F_cp3_000(a21,b3)); 
    Cp3_i00 = sign(F_cp3_i00(a21,b3)); 
    Cp3_0i0 = sign(F_cp3_0i0(a21,b3)); 
    Cp3_ii0 = sign(F_cp3_ii0(a21,b3));
    Cp3_00i = sign(F_cp3_00i(a21,b3)); 
    Cp3_i0i = sign(F_cp3_i0i(a21,b3)); 
    Cp3_0ii = sign(F_cp3_0ii(a21,b3)); 
    Cp3_iii = sign(F_cp3_iii(a21,b3));
    
    Cp4_000 = sign(F_cp4_000(a21,b3)); 
    Cp4_i00 = sign(F_cp4_i00(a21,b3)); 
    Cp4_0i0 = sign(F_cp4_0i0(a21,b3)); 
    Cp4_ii0 = sign(F_cp4_ii0(a21,b3));
    Cp4_00i = sign(F_cp4_00i(a21,b3)); 
    Cp4_i0i = sign(F_cp4_i0i(a21,b3)); 
    Cp4_0ii = sign(F_cp4_0ii(a21,b3)); 
    Cp4_iii = sign(F_cp4_iii(a21,b3));
    
    Cp5_000 = sign(F_cp5_000(a21,b3)); 
    Cp5_i00 = sign(F_cp5_i00(a21,b3)); 
    Cp5_0i0 = sign(F_cp5_0i0(a21,b3)); 
    Cp5_ii0 = sign(F_cp5_ii0(a21,b3));
    Cp5_00i = sign(F_cp5_00i(a21,b3)); 
    Cp5_i0i = sign(F_cp5_i0i(a21,b3)); 
    Cp5_0ii = sign(F_cp5_0ii(a21,b3)); 
    Cp5_iii = sign(F_cp5_iii(a21,b3));
    
    Cp6_000 = sign(F_cp6_000(a21,b3)); 
    Cp6_i00 = sign(F_cp6_i00(a21,b3)); 
    Cp6_0i0 = sign(F_cp6_0i0(a21,b3)); 
    Cp6_ii0 = sign(F_cp6_ii0(a21,b3));
    Cp6_00i = sign(F_cp6_00i(a21,b3)); 
    Cp6_i0i = sign(F_cp6_i0i(a21,b3)); 
    Cp6_0ii = sign(F_cp6_0ii(a21,b3)); 
    Cp6_iii = sign(F_cp6_iii(a21,b3));
    
V000 = (1-sign(Cp2_000))/2 + (1-sign(Cp2_000).*sign(Cp3_000))/2 + (1-sign(Cp3_000).*sign(Cp4_000))/2 + (1-sign(Cp4_000).*sign(Cp5_000))/2 + (1-sign(Cp5_000).*sign(Cp6_000))/2;
Vi00 = (1-sign(Cp2_i00))/2 + (1-sign(Cp2_i00).*sign(Cp3_i00))/2 + (1-sign(Cp3_i00).*sign(Cp4_i00))/2 + (1-sign(Cp4_i00).*sign(Cp5_i00))/2 + (1-sign(Cp5_i00).*sign(Cp6_i00))/2; 
V0i0 = (1-sign(Cp2_0i0))/2 + (1-sign(Cp2_0i0).*sign(Cp3_0i0))/2 + (1-sign(Cp3_0i0).*sign(Cp4_0i0))/2 + (1-sign(Cp4_0i0).*sign(Cp5_0i0))/2 + (1-sign(Cp5_0i0).*sign(Cp6_0i0))/2; 
V00i = (1-sign(Cp2_00i))/2 + (1-sign(Cp2_00i).*sign(Cp3_00i))/2 + (1-sign(Cp3_00i).*sign(Cp4_00i))/2 + (1-sign(Cp4_00i).*sign(Cp5_00i))/2 + (1-sign(Cp5_00i).*sign(Cp6_00i))/2; 
Vii0 = (1-sign(Cp2_ii0))/2 + (1-sign(Cp2_ii0).*sign(Cp3_ii0))/2 + (1-sign(Cp3_ii0).*sign(Cp4_ii0))/2 + (1-sign(Cp4_ii0).*sign(Cp5_ii0))/2 + (1-sign(Cp5_ii0).*sign(Cp6_ii0))/2; 
Vi0i = (1-sign(Cp2_i0i))/2 + (1-sign(Cp2_i0i).*sign(Cp3_i0i))/2 + (1-sign(Cp3_i0i).*sign(Cp4_i0i))/2 + (1-sign(Cp4_i0i).*sign(Cp5_i0i))/2 + (1-sign(Cp5_i0i).*sign(Cp6_i0i))/2; 
V0ii = (1-sign(Cp2_0ii))/2 + (1-sign(Cp2_0ii).*sign(Cp3_0ii))/2 + (1-sign(Cp3_0ii).*sign(Cp4_0ii))/2 + (1-sign(Cp4_0ii).*sign(Cp5_0ii))/2 + (1-sign(Cp5_0ii).*sign(Cp6_0ii))/2; 
Viii = (1-sign(Cp2_iii))/2 + (1-sign(Cp2_iii).*sign(Cp3_iii))/2 + (1-sign(Cp3_iii).*sign(Cp4_iii))/2 + (1-sign(Cp4_iii).*sign(Cp5_iii))/2 + (1-sign(Cp5_iii).*sign(Cp6_iii))/2; 

NpRoots = (V000 - Vi00 - V0i0 - V00i + Vii0 + Vi0i + V0ii - Viii)/4; 

Signs = [Cp2_000, Cp2_i00, Cp2_0i0, Cp2_ii0, Cp2_00i, Cp2_i0i, Cp2_0ii, Cp2_iii, ... 
         Cp3_000, Cp3_i00, Cp3_0i0, Cp3_ii0, Cp3_00i, Cp3_i0i, Cp3_0ii, Cp3_iii, ...
         Cp4_000, Cp4_i00, Cp4_0i0, Cp4_ii0, Cp4_00i, Cp4_i0i, Cp4_0ii, Cp4_iii, ...
         Cp5_000, Cp5_i00, Cp5_0i0, Cp5_ii0, Cp5_00i, Cp5_i0i, Cp5_0ii, Cp5_iii, ...
         Cp6_000, Cp6_i00, Cp6_0i0, Cp6_ii0, Cp6_00i, Cp6_i0i, Cp6_0ii, Cp6_iii, NpRoots];

unique_Signs = unique(Signs,'rows');
[rows,~] = find(unique_Signs(:,end)>0);
unique_nonzero_Signs = sortrows(unique_Signs(rows,1:end-1))

NpRoots_fix = NpRoots; 
NpRoots_fix(find(NpRoots_fix<0.51))=0;
NpRoots_fix(find(NpRoots_fix>0.51 & NpRoots_fix<1.49))=1;
NpRoots_fix(find(NpRoots_fix>1.51 & NpRoots_fix<2.51))=2;

NPROOTS = transpose(reshape(NpRoots_fix,sqrt(sample),sqrt(sample))); 

figure(1)
imagesc(b3,a21,NPROOTS);
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('b_{3}', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Number of feasible roots', 'FontSize', 30); 
annotation('textbox', [0.6, 0.25, 0.1, 0.1], 'String', ["Orange region: #of feasible roots = 2", "White region: #of feasible roots = 0"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

figure(2)
Xf=F_cp6_0i0(a21,b3);
Xf=(Xf>0 ==1) & (Xf<0 == 0);
XF = transpose(reshape(Xf,sqrt(sample),sqrt(sample))); 
imagesc(b3,a21,XF); colorbar off
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('b_{3}', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Sign of v_0(0,\infty,0)', 'FontSize', 30); 
annotation('textbox', [0.6, 0.25, 0.1, 0.1], 'String', ["Orange region: positive sign", "White region: negative sign"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

%%
figure(3)
Xf=F_cp6_000(a21,b3);
Xf=(Xf<0 ==1) & (Xf>0 == 0);
XF = transpose(reshape(Xf,sqrt(sample),sqrt(sample))); 
imagesc(b3,a21,XF); colorbar off
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('b_{3}', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Sign of v_0(0,0,0)', 'FontSize', 30); 
annotation('textbox', [0.6, 0.25, 0.1, 0.1], 'String', ["Orange region: negative sign", "White region: positive sign"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

%% Test via PHClab


    a11 =  2; 
    a12 = -1.5; 
    a13 = -1.5; 
    a22 = 2; 
    a23 = -1.5; 
    a31 = -1.5; 
    a32 = -1; 
    a33 = 1; 
    r1 = 1.5; 
    r2 = -1.5; 
    r3 = -1.5; 
    b1 = 1; 
    b2 = -1; 
    
    a21 = 4; 
    b3 = 2.5;

set_phcpath('/Users/mohammadaladwani/Documents/MATLAB/PHClab1.0.4/phc'); 

MM = [r1,  0, 0, 0; ...
      a11, 1, 0, 0; ...
      a12, 0, 1, 0; ...
      a13, 0, 0, 1; ...
      b1,  0, 1, 1; ...
      0,   0, 0, 0; ...
      r2,  0, 0, 0; ...
      a21, 1, 0, 0; ...
      a22, 0, 1, 0; ...
      a23, 0, 0, 1; ...
      b2 , 1, 0, 1; ...
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