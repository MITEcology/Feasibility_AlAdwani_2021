clear all; close all; clc; 

%% Part 1

syms r1 r2 a11 a12 a21 a22 h N1 N2 f1 f2 x y

F1=@(n1,n2) r1 + a11*n1 + a12*n1*n2 + r1*h*n1^2 + a11*h*n1^3; 
F2=@(n1,n2) r2 + a22*n2 + (a21+r2*h)*n1^2 + a22*h*n1^2*n2; 

ResN1=@(N1,N2) det([a11*h                , r1*h                  , a11 + a12*N2          , r1           , N1*f1   ; ...
                    0                    , a11*h                 , r1*h                  , a11 + a12*N2 , f1      ; ...
                   a21 + r2*h + a22*h*N2 , 0                     , r2 + a22*N2           , 0            , N1^2*f2 ; ...
                    0                    , a21 + r2*h + a22*h*N2 , 0                     , r2 + a22*N2  , N1*f2   ; ...
                    0                    , 0                     , a21 + r2*h + a22*h*N2 , 0            , f2]);

ResN2=@(N1,N2) det([a12*N1         ,f1; ...
                    a22+a22*h*N1^2 ,f2]); 
                
[cResN1, tResN1] = coeffs(ResN1(N1,N2),[f1 f2]);
T21 = cResN1(tResN1==f1); 
T22 = cResN1(tResN1==f2); 

[cResN2, tResN2] = coeffs(ResN2(N1,N2),[f1 f2]);
T11 = cResN2(tResN2==f1); 
T12 = cResN2(tResN2==f2); 

fX = matlabFunction(expand(T21*F1(N1,N2)+T22*F2(N1,N2))); 
X =@(N2) fX(N2,a11,a12,a21,a22,h,r1,r2);

fY = matlabFunction(expand(T11*F1(N1,N2)+T12*F2(N1,N2))); 
Y = @(N1) fY(N1,a11,a12,a21,a22,h,r1,r2);

H2 = 1/Y(1/x);
L2 = taylor(H2, x, 'Order', 18);
LL2 = matlabFunction(L2); 
Resultant1 =@(x) LL2(a11,a12,a21,a22,h,r1,r2,x);

H1 = 1/X(1/y);
L1 = taylor(H1, y, 'Order', 18);
LL1 = matlabFunction(L1); 
Resultant2 =@(y) LL1(a11,a12,a21,a22,h,r1,r2,y);

t11 = matlabFunction(T11); 
fT11 =@(N1) t11(N1,a22,h);
t12 = matlabFunction(T12); 
fT12 =@(N1) t12(N1,a12);
t21 = matlabFunction(T21);
fT21 =@(N1,N2) t21(N1,N2,a11,a12,a21,a22,h,r1,r2);
t22 = matlabFunction(T22);
fT22 =@(N1,N2) t22(N1,N2,a11,a12,a21,a22,h,r1,r2);

T =@(N1,N2) det([fT11(N1) fT12(N1); fT21(N1,N2) fT22(N1,N2)]); 

Jacob = det(jacobian([F1(N1,N2), F2(N1,N2)], [N1, N2]));
fJ = matlabFunction(Jacob);
J =@(N1,N2) fJ(N1,N2,a11,a12,a21,a22,h,r1,r2);

[cTJ, tTJ] = coeffs(T(x,y)*J(x,y),[x y]);
[cx, tx] = coeffs(Resultant1(x),x);
[cy, ty] = coeffs(Resultant2(y),y);

Sigma_00 = 0; Sigma_10 = 0; Sigma_20 = 0; Sigma_30 = 0; Sigma_40 = 0; Sigma_50 = 0; Sigma_60 = 0; Sigma_70 = 0; Sigma_80 = 0; Sigma_90 = 0; 
Sigma_01 = 0; Sigma_11 = 0; Sigma_21 = 0; Sigma_31 = 0; Sigma_41 = 0; Sigma_51 = 0; Sigma_61 = 0; Sigma_71 = 0; Sigma_81 = 0; Sigma_91 = 0;
Sigma_02 = 0; Sigma_12 = 0; Sigma_22 = 0; Sigma_32 = 0; Sigma_42 = 0; Sigma_52 = 0; Sigma_62 = 0; Sigma_72 = 0; Sigma_82 = 0; Sigma_92 = 0;
Sigma_03 = 0; Sigma_13 = 0; Sigma_23 = 0; Sigma_33 = 0; Sigma_43 = 0; Sigma_53 = 0; Sigma_63 = 0; Sigma_73 = 0; Sigma_83 = 0; Sigma_93 = 0;
Sigma_04 = 0; Sigma_14 = 0; Sigma_24 = 0; Sigma_34 = 0; Sigma_44 = 0; Sigma_54 = 0; Sigma_64 = 0; Sigma_74 = 0; Sigma_84 = 0; Sigma_94 = 0;
Sigma_05 = 0; Sigma_15 = 0; Sigma_25 = 0; Sigma_35 = 0; Sigma_45 = 0; Sigma_55 = 0; Sigma_65 = 0; Sigma_75 = 0; Sigma_85 = 0; Sigma_95 = 0;
Sigma_06 = 0; Sigma_16 = 0; Sigma_26 = 0; Sigma_36 = 0; Sigma_46 = 0; Sigma_56 = 0; Sigma_66 = 0; Sigma_76 = 0; Sigma_86 = 0; Sigma_96 = 0;
Sigma_07 = 0; Sigma_17 = 0; Sigma_27 = 0; Sigma_37 = 0; Sigma_47 = 0; Sigma_57 = 0; Sigma_67 = 0; Sigma_77 = 0; Sigma_87 = 0; Sigma_97 = 0;
Sigma_08 = 0; Sigma_18 = 0; Sigma_28 = 0; Sigma_38 = 0; Sigma_48 = 0; Sigma_58 = 0; Sigma_68 = 0; Sigma_78 = 0; Sigma_88 = 0; Sigma_98 = 0;
Sigma_09 = 0; Sigma_19 = 0; Sigma_29 = 0; Sigma_39 = 0; Sigma_49 = 0; Sigma_59 = 0; Sigma_69 = 0; Sigma_79 = 0; Sigma_89 = 0; Sigma_99 = 0;


for i = 1:length(tTJ)
    
    powx = polynomialDegree(tTJ(i),x);
    powy = polynomialDegree(tTJ(i),y);
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_00 = Sigma_00 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_10 = Sigma_10 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+0)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_20 = Sigma_20 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_30 = Sigma_30 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_40 = Sigma_40 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_50 = Sigma_50 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+0)); 
    end

    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_60 = Sigma_60 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_70 = Sigma_70 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_80 = Sigma_80 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+0)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+0)))==1
        Sigma_90 = Sigma_90 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+0)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_01 = Sigma_01 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_11 = Sigma_11 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+1)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_21 = Sigma_21 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_31 = Sigma_31 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_41 = Sigma_41 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_51 = Sigma_51 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_61 = Sigma_61 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_71 = Sigma_71 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_81 = Sigma_81 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+1)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+1)))==1
        Sigma_91 = Sigma_91 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+1)); 
    end
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_02 = Sigma_02 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_12 = Sigma_12 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+2)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_22 = Sigma_22 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_32 = Sigma_32 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_42 = Sigma_42 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_52 = Sigma_52 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_62 = Sigma_62 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_72 = Sigma_72 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_82 = Sigma_82 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+2)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+2)))==1
        Sigma_92 = Sigma_92 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+2)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_03 = Sigma_03 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_13 = Sigma_13 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+3)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_23 = Sigma_23 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_33 = Sigma_33 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_43 = Sigma_43 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_53 = Sigma_53 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_63 = Sigma_63 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_73 = Sigma_73 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_83 = Sigma_83 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+3)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+3)))==1
        Sigma_93 = Sigma_93 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+3)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_04 = Sigma_04 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_14 = Sigma_14 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+4)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_24 = Sigma_24 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_34 = Sigma_34 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_44 = Sigma_44 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_54 = Sigma_54 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_64 = Sigma_64 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_74 = Sigma_74 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_84 = Sigma_84 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+4)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+4)))==1
        Sigma_94 = Sigma_94 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+4)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_05 = Sigma_05 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_15 = Sigma_15 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+5)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_25 = Sigma_25 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_35 = Sigma_35 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_45 = Sigma_45 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_55 = Sigma_55 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_65 = Sigma_65 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_75 = Sigma_75 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_85 = Sigma_85 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+5)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+5)))==1
        Sigma_95 = Sigma_95 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+5)); 
    end
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_06 = Sigma_06 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_16 = Sigma_16 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+6)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_26 = Sigma_26 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_36 = Sigma_36 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_46 = Sigma_46 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_56 = Sigma_56 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_66 = Sigma_66 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_76 = Sigma_76 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_86 = Sigma_86 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+6)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+6)))==1
        Sigma_96 = Sigma_96 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+6)); 
    end  
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_07 = Sigma_07 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_17 = Sigma_17 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+7)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_27 = Sigma_27 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_37 = Sigma_37 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_47 = Sigma_47 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_57 = Sigma_57 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_67 = Sigma_67 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_77 = Sigma_77 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_87 = Sigma_87 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+7)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+7)))==1
        Sigma_97 = Sigma_97 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+7)); 
    end 
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_08 = Sigma_08 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_18 = Sigma_18 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+8)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_28 = Sigma_28 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_38 = Sigma_38 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_48 = Sigma_48 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_58 = Sigma_58 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_68 = Sigma_68 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_78 = Sigma_78 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_88 = Sigma_88 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+8)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+8)))==1
        Sigma_98 = Sigma_98 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+8)); 
    end 
    
    %
    
    if max(ismember(tx,x^(powx+1+0)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_09 = Sigma_09 + cTJ(i)*cx(tx==x^(powx+1+0))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+1)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_19 = Sigma_19 + cTJ(i)*cx(tx==x^(powx+1+1))*cy(ty==y^(powy+1+9)); 
    end

    if max(ismember(tx,x^(powx+1+2)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_29 = Sigma_29 + cTJ(i)*cx(tx==x^(powx+1+2))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+3)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_39 = Sigma_39 + cTJ(i)*cx(tx==x^(powx+1+3))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+4)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_49 = Sigma_49 + cTJ(i)*cx(tx==x^(powx+1+4))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+5)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_59 = Sigma_59 + cTJ(i)*cx(tx==x^(powx+1+5))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+6)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_69 = Sigma_69 + cTJ(i)*cx(tx==x^(powx+1+6))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+7)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_79 = Sigma_79 + cTJ(i)*cx(tx==x^(powx+1+7))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+8)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_89 = Sigma_89 + cTJ(i)*cx(tx==x^(powx+1+8))*cy(ty==y^(powy+1+9)); 
    end
    
    if max(ismember(tx,x^(powx+1+9)))*max(ismember(ty,y^(powy+1+9)))==1
        Sigma_99 = Sigma_99 + cTJ(i)*cx(tx==x^(powx+1+9))*cy(ty==y^(powy+1+9)); 
    end    
    
end

%% part 2

syms N11 N12 N13 N14 N15 N21 N22 N23 N24 N25 s1 s2 ...
    c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14

c1 = 1;  c2 = 0; c3 = 0; c4 = 1; c5 = 0; c6 = 0; c7 = 1; c8 = 0; c9 = 0; 
c10= 0;  c11= 0; c12= 1; c13= 0; c14= 0; 

%m=@(N1,N2) [1; c1*N1+c2*N2; c3*N1^2+c4*N1*N2+c5*N2^2; c6*N1^3+c7*N1^2*N2+c8*N1*N2^2+c9*N2^3; c10*N1^4+c11*N1^3*N2+c12*N1^2*N2^2+c13*N1*N2^3+c14*N2^4];
m=@(N1,N2) [1; N1; N2; N1.*N2; N1.^2];

W=[m(N11,N21),m(N12,N22),m(N13,N23),m(N14,N24),m(N15,N25)]; 
D = diag([(N11-s1)*(N21-s2),(N12-s1)*(N22-s2),(N13-s1)*(N23-s2),(N14-s1)*(N24-s2),(N15-s1)*(N25-s2)]);
S=collect(W*D*transpose(W),[s1 s2 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14]); 

syms      S10 S20 S30 S40 S50 S60 S70 S80 S90 ...
      S01 S11 S21 S31 S41 S51 S61 S71 S81 S91 ...
      S02 S12 S22 S32 S42 S52 S62 S72 S82 S92 ...
      S03 S13 S23 S33 S43 S53 S63 S73 S83 S93 ...
      S04 S14 S24 S34 S44 S54 S64 S74 S84 S94 ...
      S05 S15 S25 S35 S45 S55 S65 S75 S85 S95 ...
      S06 S16 S26 S36 S46 S56 S66 S76 S86 S96 ...
      S07 S17 S27 S37 S47 S57 S67 S77 S87 S97 ...
      S08 S18 S28 S38 S48 S58 S68 S78 S88 S98 ...
      S09 S19 S29 S39 S49 S59 S69 S79 S89 S99
 
  
sums =[...
N11^1*N21^0 + N12^1*N22^0 + N13^1*N23^0 + N14^1*N24^0 + N15^1*N25^0, ... %
N11^2*N21^0 + N12^2*N22^0 + N13^2*N23^0 + N14^2*N24^0 + N15^2*N25^0, ...
N11^3*N21^0 + N12^3*N22^0 + N13^3*N23^0 + N14^3*N24^0 + N15^3*N25^0, ...
N11^4*N21^0 + N12^4*N22^0 + N13^4*N23^0 + N14^4*N24^0 + N15^4*N25^0, ...
N11^5*N21^0 + N12^5*N22^0 + N13^5*N23^0 + N14^5*N24^0 + N15^5*N25^0, ...
N11^6*N21^0 + N12^6*N22^0 + N13^6*N23^0 + N14^6*N24^0 + N15^6*N25^0, ...
N11^7*N21^0 + N12^7*N22^0 + N13^7*N23^0 + N14^7*N24^0 + N15^7*N25^0, ...
N11^8*N21^0 + N12^8*N22^0 + N13^8*N23^0 + N14^8*N24^0 + N15^8*N25^0, ...
N11^9*N21^0 + N12^9*N22^0 + N13^9*N23^0 + N14^9*N24^0 + N15^9*N25^0, ... 
N11^0*N21^1 + N12^0*N22^1 + N13^0*N23^1 + N14^0*N24^1 + N15^0*N25^1, ... % 
N11^1*N21^1 + N12^1*N22^1 + N13^1*N23^1 + N14^1*N24^1 + N15^1*N25^1, ...
N11^2*N21^1 + N12^2*N22^1 + N13^2*N23^1 + N14^2*N24^1 + N15^2*N25^1, ...
N11^3*N21^1 + N12^3*N22^1 + N13^3*N23^1 + N14^3*N24^1 + N15^3*N25^1, ...
N11^4*N21^1 + N12^4*N22^1 + N13^4*N23^1 + N14^4*N24^1 + N15^4*N25^1, ...
N11^5*N21^1 + N12^5*N22^1 + N13^5*N23^1 + N14^5*N24^1 + N15^5*N25^1, ...
N11^6*N21^1 + N12^6*N22^1 + N13^6*N23^1 + N14^6*N24^1 + N15^6*N25^1, ...
N11^7*N21^1 + N12^7*N22^1 + N13^7*N23^1 + N14^7*N24^1 + N15^7*N25^1, ...
N11^8*N21^1 + N12^8*N22^1 + N13^8*N23^1 + N14^8*N24^1 + N15^8*N25^1, ...
N11^9*N21^1 + N12^9*N22^1 + N13^9*N23^1 + N14^9*N24^1 + N15^9*N25^1, ...
N11^0*N21^2 + N12^0*N22^2 + N13^0*N23^2 + N14^0*N24^2 + N15^0*N25^2, ... % 
N11^1*N21^2 + N12^1*N22^2 + N13^1*N23^2 + N14^1*N24^2 + N15^1*N25^2, ...
N11^2*N21^2 + N12^2*N22^2 + N13^2*N23^2 + N14^2*N24^2 + N15^2*N25^2, ...
N11^3*N21^2 + N12^3*N22^2 + N13^3*N23^2 + N14^3*N24^2 + N15^3*N25^2, ...
N11^4*N21^2 + N12^4*N22^2 + N13^4*N23^2 + N14^4*N24^2 + N15^4*N25^2, ...
N11^5*N21^2 + N12^5*N22^2 + N13^5*N23^2 + N14^5*N24^2 + N15^5*N25^2, ...
N11^6*N21^2 + N12^6*N22^2 + N13^6*N23^2 + N14^6*N24^2 + N15^6*N25^2, ...
N11^7*N21^2 + N12^7*N22^2 + N13^7*N23^2 + N14^7*N24^2 + N15^7*N25^2, ...
N11^8*N21^2 + N12^8*N22^2 + N13^8*N23^2 + N14^8*N24^2 + N15^8*N25^2, ...
N11^9*N21^2 + N12^9*N22^2 + N13^9*N23^2 + N14^9*N24^2 + N15^9*N25^2, ...
N11^0*N21^3 + N12^0*N22^3 + N13^0*N23^3 + N14^0*N24^3 + N15^0*N25^3, ... % 
N11^1*N21^3 + N12^1*N22^3 + N13^1*N23^3 + N14^1*N24^3 + N15^1*N25^3, ...
N11^2*N21^3 + N12^2*N22^3 + N13^2*N23^3 + N14^2*N24^3 + N15^2*N25^3, ...
N11^3*N21^3 + N12^3*N22^3 + N13^3*N23^3 + N14^3*N24^3 + N15^3*N25^3, ...
N11^4*N21^3 + N12^4*N22^3 + N13^4*N23^3 + N14^4*N24^3 + N15^4*N25^3, ...
N11^5*N21^3 + N12^5*N22^3 + N13^5*N23^3 + N14^5*N24^3 + N15^5*N25^3, ...
N11^6*N21^3 + N12^6*N22^3 + N13^6*N23^3 + N14^6*N24^3 + N15^6*N25^3, ...
N11^7*N21^3 + N12^7*N22^3 + N13^7*N23^3 + N14^7*N24^3 + N15^7*N25^3, ...
N11^8*N21^3 + N12^8*N22^3 + N13^8*N23^3 + N14^8*N24^3 + N15^8*N25^3, ...
N11^9*N21^3 + N12^9*N22^3 + N13^9*N23^3 + N14^9*N24^3 + N15^9*N25^3, ...
N11^0*N21^4 + N12^0*N22^4 + N13^0*N23^4 + N14^0*N24^4 + N15^0*N25^4, ... % 
N11^1*N21^4 + N12^1*N22^4 + N13^1*N23^4 + N14^1*N24^4 + N15^1*N25^4, ...
N11^2*N21^4 + N12^2*N22^4 + N13^2*N23^4 + N14^2*N24^4 + N15^2*N25^4, ...
N11^3*N21^4 + N12^3*N22^4 + N13^3*N23^4 + N14^3*N24^4 + N15^3*N25^4, ...
N11^4*N21^4 + N12^4*N22^4 + N13^4*N23^4 + N14^4*N24^4 + N15^4*N25^4, ...
N11^5*N21^4 + N12^5*N22^4 + N13^5*N23^4 + N14^5*N24^4 + N15^5*N25^4, ...
N11^6*N21^4 + N12^6*N22^4 + N13^6*N23^4 + N14^6*N24^4 + N15^6*N25^4, ...
N11^7*N21^4 + N12^7*N22^4 + N13^7*N23^4 + N14^7*N24^4 + N15^7*N25^4, ...
N11^8*N21^4 + N12^8*N22^4 + N13^8*N23^4 + N14^8*N24^4 + N15^8*N25^4, ...
N11^9*N21^4 + N12^9*N22^4 + N13^9*N23^4 + N14^9*N24^4 + N15^9*N25^4, ...
N11^0*N21^5 + N12^0*N22^5 + N13^0*N23^5 + N14^0*N24^5 + N15^0*N25^5, ... % 
N11^1*N21^5 + N12^1*N22^5 + N13^1*N23^5 + N14^1*N24^5 + N15^1*N25^5, ...
N11^2*N21^5 + N12^2*N22^5 + N13^2*N23^5 + N14^2*N24^5 + N15^2*N25^5, ...
N11^3*N21^5 + N12^3*N22^5 + N13^3*N23^5 + N14^3*N24^5 + N15^3*N25^5, ...
N11^4*N21^5 + N12^4*N22^5 + N13^4*N23^5 + N14^4*N24^5 + N15^4*N25^5, ...
N11^5*N21^5 + N12^5*N22^5 + N13^5*N23^5 + N14^5*N24^5 + N15^5*N25^5, ...
N11^6*N21^5 + N12^6*N22^5 + N13^6*N23^5 + N14^6*N24^5 + N15^6*N25^5, ...
N11^7*N21^5 + N12^7*N22^5 + N13^7*N23^5 + N14^7*N24^5 + N15^7*N25^5, ...
N11^8*N21^5 + N12^8*N22^5 + N13^8*N23^5 + N14^8*N24^5 + N15^8*N25^5, ...
N11^9*N21^5 + N12^9*N22^5 + N13^9*N23^5 + N14^9*N24^5 + N15^9*N25^5, ...
N11^0*N21^6 + N12^0*N22^6 + N13^0*N23^6 + N14^0*N24^6 + N15^0*N25^6, ... % 
N11^1*N21^6 + N12^1*N22^6 + N13^1*N23^6 + N14^1*N24^6 + N15^1*N25^6, ...
N11^2*N21^6 + N12^2*N22^6 + N13^2*N23^6 + N14^2*N24^6 + N15^2*N25^6, ...
N11^3*N21^6 + N12^3*N22^6 + N13^3*N23^6 + N14^3*N24^6 + N15^3*N25^6, ...
N11^4*N21^6 + N12^4*N22^6 + N13^4*N23^6 + N14^4*N24^6 + N15^4*N25^6, ...
N11^5*N21^6 + N12^5*N22^6 + N13^5*N23^6 + N14^5*N24^6 + N15^5*N25^6, ...
N11^6*N21^6 + N12^6*N22^6 + N13^6*N23^6 + N14^6*N24^6 + N15^6*N25^6, ...
N11^7*N21^6 + N12^7*N22^6 + N13^7*N23^6 + N14^7*N24^6 + N15^7*N25^6, ...
N11^8*N21^6 + N12^8*N22^6 + N13^8*N23^6 + N14^8*N24^6 + N15^8*N25^6, ...
N11^9*N21^6 + N12^9*N22^6 + N13^9*N23^6 + N14^9*N24^6 + N15^9*N25^6, ...
N11^0*N21^7 + N12^0*N22^7 + N13^0*N23^7 + N14^0*N24^7 + N15^0*N25^7, ... % 
N11^1*N21^7 + N12^1*N22^7 + N13^1*N23^7 + N14^1*N24^7 + N15^1*N25^7, ...
N11^2*N21^7 + N12^2*N22^7 + N13^2*N23^7 + N14^2*N24^7 + N15^2*N25^7, ...
N11^3*N21^7 + N12^3*N22^7 + N13^3*N23^7 + N14^3*N24^7 + N15^3*N25^7, ...
N11^4*N21^7 + N12^4*N22^7 + N13^4*N23^7 + N14^4*N24^7 + N15^4*N25^7, ...
N11^5*N21^7 + N12^5*N22^7 + N13^5*N23^7 + N14^5*N24^7 + N15^5*N25^7, ...
N11^6*N21^7 + N12^6*N22^7 + N13^6*N23^7 + N14^6*N24^7 + N15^6*N25^7, ...
N11^7*N21^7 + N12^7*N22^7 + N13^7*N23^7 + N14^7*N24^7 + N15^7*N25^7, ...
N11^8*N21^7 + N12^8*N22^7 + N13^8*N23^7 + N14^8*N24^7 + N15^8*N25^7, ...
N11^9*N21^7 + N12^9*N22^7 + N13^9*N23^7 + N14^9*N24^7 + N15^9*N25^7, ...
N11^0*N21^8 + N12^0*N22^8 + N13^0*N23^8 + N14^0*N24^8 + N15^0*N25^8, ... % 
N11^1*N21^8 + N12^1*N22^8 + N13^1*N23^8 + N14^1*N24^8 + N15^1*N25^8, ...
N11^2*N21^8 + N12^2*N22^8 + N13^2*N23^8 + N14^2*N24^8 + N15^2*N25^8, ...
N11^3*N21^8 + N12^3*N22^8 + N13^3*N23^8 + N14^3*N24^8 + N15^3*N25^8, ...
N11^4*N21^8 + N12^4*N22^8 + N13^4*N23^8 + N14^4*N24^8 + N15^4*N25^8, ...
N11^5*N21^8 + N12^5*N22^8 + N13^5*N23^8 + N14^5*N24^8 + N15^5*N25^8, ...
N11^6*N21^8 + N12^6*N22^8 + N13^6*N23^8 + N14^6*N24^8 + N15^6*N25^8, ...
N11^7*N21^8 + N12^7*N22^8 + N13^7*N23^8 + N14^7*N24^8 + N15^7*N25^8, ...
N11^8*N21^8 + N12^8*N22^8 + N13^8*N23^8 + N14^8*N24^8 + N15^8*N25^8, ...
N11^9*N21^8 + N12^9*N22^8 + N13^9*N23^8 + N14^9*N24^8 + N15^9*N25^8, ...
N11^0*N21^9 + N12^0*N22^9 + N13^0*N23^9 + N14^0*N24^9 + N15^0*N25^9, ... % 
N11^1*N21^9 + N12^1*N22^9 + N13^1*N23^9 + N14^1*N24^9 + N15^1*N25^9, ...
N11^2*N21^9 + N12^2*N22^9 + N13^2*N23^9 + N14^2*N24^9 + N15^2*N25^9, ...
N11^3*N21^9 + N12^3*N22^9 + N13^3*N23^9 + N14^3*N24^9 + N15^3*N25^9, ...
N11^4*N21^9 + N12^4*N22^9 + N13^4*N23^9 + N14^4*N24^9 + N15^4*N25^9, ...
N11^5*N21^9 + N12^5*N22^9 + N13^5*N23^9 + N14^5*N24^9 + N15^5*N25^9, ...
N11^6*N21^9 + N12^6*N22^9 + N13^6*N23^9 + N14^6*N24^9 + N15^6*N25^9, ...
N11^7*N21^9 + N12^7*N22^9 + N13^7*N23^9 + N14^7*N24^9 + N15^7*N25^9, ...
N11^8*N21^9 + N12^8*N22^9 + N13^8*N23^9 + N14^8*N24^9 + N15^8*N25^9, ...
N11^9*N21^9 + N12^9*N22^9 + N13^9*N23^9 + N14^9*N24^9 + N15^9*N25^9];

notions = [...
          S10 S20 S30 S40 S50 S60 S70 S80 S90 ...
      S01 S11 S21 S31 S41 S51 S61 S71 S81 S91 ...
      S02 S12 S22 S32 S42 S52 S62 S72 S82 S92 ...
      S03 S13 S23 S33 S43 S53 S63 S73 S83 S93 ...
      S04 S14 S24 S34 S44 S54 S64 S74 S84 S94 ...
      S05 S15 S25 S35 S45 S55 S65 S75 S85 S95 ...
      S06 S16 S26 S36 S46 S56 S66 S76 S86 S96 ...
      S07 S17 S27 S37 S47 S57 S67 S77 S87 S97 ...
      S08 S18 S28 S38 S48 S58 S68 S78 S88 S98 ...
      S09 S19 S29 S39 S49 S59 S69 S79 S89 S99];
        
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

[cp1,tp1] = coeffs(P(1),[s1 s2]); 
[cp2,tp2] = coeffs(P(2),[s1 s2]); 
[cp3,tp3] = coeffs(P(3),[s1 s2]);
[cp4,tp4] = coeffs(P(4),[s1 s2]); 
[cp5,tp5] = coeffs(P(5),[s1 s2]);
[cp6,tp6] = coeffs(P(6),[s1 s2]);

cp1_00 = 1;
cp1_i0 = 1; 
cp1_0i = 1;
cp1_ii = 1;

cp2_00 = cp2(tp2==1);
cp2_i0 = cp2(tp2==s1); 
cp2_0i = cp2(tp2==s2);
cp2_ii = cp2(tp2==s1*s2); 

cp3_00 = cp3(tp3==1);
cp3_i0 = cp3(tp3==s1^2); 
cp3_0i = cp3(tp3==s2^2);
cp3_ii = cp3(tp3==s1^2*s2^2); 

cp4_00 = cp4(tp4==1);
cp4_i0 = cp4(tp4==s1^3); 
cp4_0i = cp4(tp4==s2^3);
cp4_ii = cp4(tp4==s1^3*s2^3); 

cp5_00 = cp5(tp5==1);
cp5_i0 = cp5(tp5==s1^4); 
cp5_0i = cp5(tp5==s2^4);
cp5_ii = cp5(tp5==s1^4*s2^4); 

cp6_00 = cp6(tp6==1);
cp6_i0 = cp6(tp6==s1^5); 
cp6_0i = cp6(tp6==s2^5);
cp6_ii = cp6(tp6==s1^5*s2^5); 

                        S10=simplify(Sigma_10); S20=simplify(Sigma_20); S30=simplify(Sigma_30); S40=simplify(Sigma_40); S50=simplify(Sigma_50); S60=simplify(Sigma_60); S70=simplify(Sigma_70); S80=simplify(Sigma_80); S90=simplify(Sigma_90);
S01=simplify(Sigma_01); S11=simplify(Sigma_11); S21=simplify(Sigma_21); S31=simplify(Sigma_31); S41=simplify(Sigma_41); S51=simplify(Sigma_51); S61=simplify(Sigma_61); S71=simplify(Sigma_71); S81=simplify(Sigma_81); S91=simplify(Sigma_91); 
S02=simplify(Sigma_02); S12=simplify(Sigma_12); S22=simplify(Sigma_22); S32=simplify(Sigma_32); S42=simplify(Sigma_42); S52=simplify(Sigma_52); S62=simplify(Sigma_62); S72=simplify(Sigma_72); S82=simplify(Sigma_82); S92=simplify(Sigma_92); 
S03=simplify(Sigma_03); S13=simplify(Sigma_13); S23=simplify(Sigma_23); S33=simplify(Sigma_33); S43=simplify(Sigma_43); S53=simplify(Sigma_53); S63=simplify(Sigma_63); S73=simplify(Sigma_73); S83=simplify(Sigma_83); S93=simplify(Sigma_93); 
S04=simplify(Sigma_04); S14=simplify(Sigma_14); S24=simplify(Sigma_24); S34=simplify(Sigma_34); S44=simplify(Sigma_44); S54=simplify(Sigma_54); S64=simplify(Sigma_64); S74=simplify(Sigma_74); S84=simplify(Sigma_84); S94=simplify(Sigma_94); 
S05=simplify(Sigma_05); S15=simplify(Sigma_15); S25=simplify(Sigma_25); S35=simplify(Sigma_35); S45=simplify(Sigma_45); S55=simplify(Sigma_55); S65=simplify(Sigma_65); S75=simplify(Sigma_75); S85=simplify(Sigma_85); S95=simplify(Sigma_95); 
S06=simplify(Sigma_06); S16=simplify(Sigma_16); S26=simplify(Sigma_26); S36=simplify(Sigma_36); S46=simplify(Sigma_46); S56=simplify(Sigma_56); S66=simplify(Sigma_66); S76=simplify(Sigma_76); S86=simplify(Sigma_86); S96=simplify(Sigma_96); 
S07=simplify(Sigma_07); S17=simplify(Sigma_17); S27=simplify(Sigma_27); S37=simplify(Sigma_37); S47=simplify(Sigma_47); S57=simplify(Sigma_57); S67=simplify(Sigma_67); S77=simplify(Sigma_77); S87=simplify(Sigma_87); S97=simplify(Sigma_97); 
S08=simplify(Sigma_08); S18=simplify(Sigma_18); S28=simplify(Sigma_28); S38=simplify(Sigma_38); S48=simplify(Sigma_48); S58=simplify(Sigma_58); S68=simplify(Sigma_68); S78=simplify(Sigma_78); S88=simplify(Sigma_88); S98=simplify(Sigma_98); 
S09=simplify(Sigma_09); S19=simplify(Sigma_19); S29=simplify(Sigma_29); S39=simplify(Sigma_39); S49=simplify(Sigma_49); S59=simplify(Sigma_59); S69=simplify(Sigma_69); S79=simplify(Sigma_79); S89=simplify(Sigma_89); S99=simplify(Sigma_99); 

F_cp2_00 = matlabFunction(subs(cp2_00));
F_cp2_i0 = matlabFunction(subs(cp2_i0));
F_cp2_0i = matlabFunction(subs(cp2_0i));
F_cp2_ii = matlabFunction(subs(cp2_ii));

F_cp3_00 = matlabFunction(subs(cp3_00));
F_cp3_i0 = matlabFunction(subs(cp3_i0));
F_cp3_0i = matlabFunction(subs(cp3_0i));
F_cp3_ii = matlabFunction(subs(cp3_ii));

F_cp4_00 = matlabFunction(subs(cp4_00));
F_cp4_i0 = matlabFunction(subs(cp4_i0));
F_cp4_0i = matlabFunction(subs(cp4_0i));
F_cp4_ii = matlabFunction(subs(cp4_ii));

F_cp5_00 = matlabFunction(subs(cp5_00));
F_cp5_i0 = matlabFunction(subs(cp5_i0));
F_cp5_0i = matlabFunction(subs(cp5_0i));
F_cp5_ii = matlabFunction(subs(cp5_ii));

F_cp6_00 = matlabFunction(subs(cp6_00));
F_cp6_i0 = matlabFunction(subs(cp6_i0));
F_cp6_0i = matlabFunction(subs(cp6_0i));
F_cp6_ii = matlabFunction(subs(cp6_ii));

%% Example

sample = 2^16; 

    a11 = 1*ones(sample,1); 
    a12 = -1.5*ones(sample,1);
    a21 = linspace(-6,-1,sqrt(sample))'; 
    a22 = 1*ones(sample,1); 
     r1 =  0.5*ones(sample,1); 
     r2 = -1.5*ones(sample,1); 
     h = linspace(0.5,4,sqrt(sample))'; % 2 equilibria
       
     [A,H] = meshgrid(a21,h); 
     a21 = reshape(A,sample,1); 
     h = reshape(H,sample,1); 
     
    Cp2_00 = sign(F_cp2_00(a11,a12,a21,a22,h,r1,r2)); 
    Cp2_i0 = sign(F_cp2_i0(a11,a12,a21,a22,h,r1,r2));
    Cp2_0i = sign(F_cp2_0i(a11,a12,a21,a22,h,r1,r2));
    Cp2_ii = sign(F_cp2_ii(a11,a12,a21,a22,h,r1,r2));  
    Cp3_00 = sign(F_cp3_00(a11,a12,a21,a22,h,r1,r2)); 
    Cp3_i0 = sign(F_cp3_i0(a11,a12,a21,a22,h,r1,r2)); 
    Cp3_0i = sign(F_cp3_0i(a11,a12,a21,a22,h,r1,r2)); 
    Cp3_ii = sign(F_cp3_ii(a11,a12,a21,a22,h,r1,r2));
    Cp4_00 = sign(F_cp4_00(a11,a12,a21,a22,h,r1,r2)); 
    Cp4_i0 = sign(F_cp4_i0(a11,a12,a21,a22,h,r1,r2)); 
    Cp4_0i = sign(F_cp4_0i(a11,a12,a21,a22,h,r1,r2)); 
    Cp4_ii = sign(F_cp4_ii(a11,a12,a21,a22,h,r1,r2)); 
    Cp5_00 = sign(F_cp5_00(a11,a12,a21,a22,h,r1,r2)); 
    Cp5_i0 = sign(F_cp5_i0(a11,a12,a21,a22,h,r1,r2)); 
    Cp5_0i = sign(F_cp5_0i(a11,a12,a21,a22,h,r1,r2)); 
    Cp5_ii = sign(F_cp5_ii(a11,a12,a21,a22,h,r1,r2)); 
    Cp6_00 = sign(F_cp6_00(a11,a12,a21,a22,h,r1,r2)); 
    Cp6_i0 = sign(F_cp6_i0(a11,a12,a21,a22,h,r1,r2)); 
    Cp6_0i = sign(F_cp6_0i(a11,a12,a21,a22,h,r1,r2)); 
    Cp6_ii = sign(F_cp6_ii(a11,a12,a21,a22,h,r1,r2)); 
     
V00 = (1-(Cp2_00))/2 + (1-(Cp2_00).*(Cp3_00))/2 + (1-(Cp3_00).*(Cp4_00))/2 + (1-(Cp4_00).*(Cp5_00))/2 + (1-(Cp5_00).*(Cp6_00))/2;
Vi0 = (1-(Cp2_i0))/2 + (1-(Cp2_i0).*(Cp3_i0))/2 + (1-(Cp3_i0).*(Cp4_i0))/2 + (1-(Cp4_i0).*(Cp5_i0))/2 + (1-(Cp5_i0).*(Cp6_i0))/2; 
V0i = (1-(Cp2_0i))/2 + (1-(Cp2_0i).*(Cp3_0i))/2 + (1-(Cp3_0i).*(Cp4_0i))/2 + (1-(Cp4_0i).*(Cp5_0i))/2 + (1-(Cp5_0i).*(Cp6_0i))/2; 
Vii = (1-(Cp2_ii))/2 + (1-(Cp2_ii).*(Cp3_ii))/2 + (1-(Cp3_ii).*(Cp4_ii))/2 + (1-(Cp4_ii).*(Cp5_ii))/2 + (1-(Cp5_ii).*(Cp6_ii))/2;

NpRoots = (V00 - Vi0 - V0i + Vii)/2; 

Signs = [Cp2_00, Cp2_i0, Cp2_0i, Cp2_ii, Cp3_00, Cp3_i0, Cp3_0i, Cp3_ii, Cp4_00, Cp4_i0, Cp4_0i, Cp4_ii, ...
         Cp5_00, Cp5_i0, Cp5_0i, Cp5_ii, Cp6_00, Cp6_i0, Cp6_0i, Cp6_ii, NpRoots]

unique_Signs = unique(Signs,'rows');
[rows,~] = find(unique_Signs(:,end)>0);
unique_nonzero_Signs = sortrows(unique_Signs(rows,1:end-1));

NPROOTS = transpose(reshape(NpRoots,sqrt(sample),sqrt(sample))); 
figure(1)
imagesc(linspace(0.5,4,sqrt(sample))',linspace(-6,-1,sqrt(sample))',NPROOTS);
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('h', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Number of feasible roots', 'FontSize', 30); 
annotation('textbox', [0.6, 0.2, 0.1, 0.1], 'String', ["Orange region: #of feasible roots = 2", "White region: #of feasible roots = 0"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

figure(2)
Xf=F_cp6_00(a11,a12,a21,a22,h,r1,r2);
Xf=(Xf<0 ==1) & (Xf>0 == 0);
XF = transpose(reshape(Xf,sqrt(sample),sqrt(sample))); 
imagesc(h,a21,XF); colorbar off
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('h', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Sign of v_0(0,0)', 'FontSize', 30); 
annotation('textbox', [0.6, 0.2, 0.1, 0.1], 'String', ["Orange region: negative sign", "White region: positive sign"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

%% Check Numerical Solution

    a11 = 1.0; 
    a12 = -1.5;
    a21 = -5.7; 
    a22 = 1; 
     r1 = 0.5; 
     r2 = -1.5; 
     h = 3; % 2 equilibria

set_phcpath('/Users/mohammadaladwani/Documents/MATLAB/PHClab1.0.4/phc'); 

MM = [r1,       0, 0; ...
      a11,      1, 0; ...
      a12,      1, 1; ...
      r1*h,     2, 0; ...
      a11*h,    3, 0; ...
      0,        0, 0; ...
      r2,       0, 0; ...
      a22,      0, 1; ...
      a21+r2*h, 2, 0; ...
      a22*h,    2, 1; ...
      0,        0, 0];    
  
make_system(MM); % shows symbolic format of the system
s = solve_system(MM); % call the blackbox solver
ns = size(s,2); % check the number of solutions
 
X = [s.x1; s.x2]
count_feasible = 0;
size(X,2);

for i = 1:size(X,2)
    
    if (real(X(1,i)) > 10^-6 && real(X(2,i)) > 10^-6 && ... 
    abs(imag(X(1,i))) < 10^-6 && abs(imag(X(2,i))) < 10^-6)

         count_feasible  = count_feasible  + 1;
    
    end

    
end

Number_of_feasible_roots = count_feasible

%% Check Feasibility Domain 

sample = 2^12; 

    a11 = 1.0; 
    a12 = -1.5;
    a22 = 1; 
     r1 = 0.5; 
     r2 = -1.5; 

    a21 = linspace(-6,-1,sqrt(sample))'; 
    h = linspace(0.5,4,sqrt(sample))'; % 2 equilibria
       
     [A,H] = meshgrid(a21,h); 
     a21 = reshape(A,sample,1); 
     h = reshape(H,sample,1); 

set_phcpath('/Users/mohammadaladwani/Documents/MATLAB/PHClab1.0.4/phc'); 

for k = 1:sample
    
    k

MM = [r1,             0, 0; ...
      a11,            1, 0; ...
      a12,            1, 1; ...
      r1*h(k),        2, 0; ...
      a11*h(k),       3, 0; ...
      0,              0, 0; ...
      r2,             0, 0; ...
      a22,            0, 1; ...
      a21(k)+r2*h(k), 2, 0; ...
      a22*h(k),       2, 1; ...
      0,              0, 0];    
  
make_system(MM); % shows symbolic format of the system
s = solve_system(MM); % call the blackbox solver
ns = size(s,2); % check the number of solutions
 
X = [s.x1; s.x2]; 
count_feasible(k) = 0;
size(X,2);

for i = 1:size(X,2)
    
    if (real(X(1,i)) > 10^-6 && real(X(2,i)) > 10^-6 && ... 
    abs(imag(X(1,i))) < 10^-6 && abs(imag(X(2,i))) < 10^-6)

         count_feasible(k)  = count_feasible(k)  + 1;
    
    end

    
end

end

count_feasible = reshape(count_feasible,sqrt(sample),sqrt(sample)); 
figure(3)
imagesc(linspace(0.5,4,sqrt(sample))',linspace(-6,-1,sqrt(sample))',count_feasible);
mycolors = [1 1 1; 0.9100 0.4100 0.1700];
colormap(mycolors);
xlabel('h', 'FontSize', 30); 
ylabel('a_{21}', 'FontSize', 30); 
title('Number of feasible roots', 'FontSize', 30); 
annotation('textbox', [0.6, 0.2, 0.1, 0.1], 'String', ["Orange region: #of feasible roots = 2", "White region: #of feasible roots = 0"], 'FontSize', 13,'FitBoxToText','on','Horizontalalignment','center')

