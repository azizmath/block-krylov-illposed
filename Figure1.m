%  graph 1 fig==1 symmetric  
clear ;clc;clf;
n=200;p=5;noiselevel=10^-2;
[A] = deriv2(n,1); [AA,b,x] = phillips(n);A=-A;
%[A,b,x] = foxgood(n); 
%A[A] = gravity(n,1,0,1,0.25); [AA,b,x] = phillips(n);
%[A,b,x] = phillips(n); 
%[A,b,x] = shaw(n); 
%[A,b,x] = heat(n,1);
%A = gallery('lotkin',n);[At b x] = shaw(n);
t=linspace(-6,6,n);
y=(1/2)*cos((1/3)*t)+(1/4);
y=y';
X_true=zeros(n,p);
X_true(:,1)=x;
for j=2:p
  X_true(:,j)=X_true(:,j-1)+(1/2)*y;
end
Btrue=A*X_true;
E=randn(n,p);
E=E/norm(E,'fro');
E=noiselevel*norm(Btrue,'fro')*E;
B=Btrue+E;
norm_error=norm(E,'fro');
m=40; fig=1;
%%%%% Gravity and Phillips fig==1
if fig==1
    [Q_m,Q_mm,T_mm,rhs,C,normB,breaks]= Block_Lanczos_tridiagonalization(A,B,m,p,1);
 else
    %%%% Nonsymmetric part 
    [Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks]= Block_Golub_Kahan_bidiagonalization(A,B,m,p,1);
 end