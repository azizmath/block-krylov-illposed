clear
n=200;p=5;noiselevel=10^-2;max_iter=40;
ylimit = [1e-16 1e1];
makefigs = 0;
reorth = 1;

%[A,b,x] = baart(n); titolo = 'Baart';
[A,b,x] = heat(n,1); titolo = 'Heat';
%A = gallery('lotkin',n);[At, b, x] = shaw(n); titolo = 'Lotkin';
%[A,b,x] = wing(n,1/3,2/3); titolo = 'Wing';

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

[U_A,Sigma_A,V_A]=csvd(A);
[Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks]= Block_Golub_Kahan_bidiagonalization(A,B,max_iter,p,reorth);
if breaks, max_iter = breaks; end

[U_C,Sigma_C,V_C]=csvd(C_mm);

R_sigma = zeros(max_iter,1);
vk = [1:max_iter]';
for k = 1:max_iter
    k1=round(p*k/3);
    if k1==0
        k1=1;
    end
    xx=zeros(k1,1);
    for j=1:k1
        xx(j,1)=abs(Sigma_C(j,1)-Sigma_A(j,1))/abs(Sigma_A(j,1));
    end
   R_sigma(k,1)=max(xx);
end

figure(1)
semilogy(vk, R_sigma, 'o-')
set(gca,'fontsize',12)
xlim(vk([1 end]))
ylim(ylimit)
title(titolo);
grid
