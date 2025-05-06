% This code provides us with the results for  plots of the Fig 3
clear; clc;clf;
n=200;p=5;noiselevel=0;max_iter=11;
ylimit = [1e-16 1e1];
makefigs = 0;
reorth = 1;

%[A] = deriv2(n,1); [AA,b,x] = phillips(n);A=-A; titolo = 'Deriv2';
%[A,b,x] = foxgood(n); titolo = 'Foxgood';
[A] = gravity(n,1,0,1,0.25); [AA,b,x] = phillips(n); titolo = 'Gravity';
%[A,b,x] = phillips(n); titolo = 'Phillips';


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

[W,LambdaA]=eig(A);
lambdaA=diag(LambdaA);
[sa vs] = sort(abs(lambdaA),'descend');
W = W(:,vs);
lambdaA= lambdaA(vs);
%%% Aziz starts here for figure 2 in papaer 1
 [Q_k,Q_kk,T_kk,rhs,C,normB,breaks]= Block_Lanczos_tridiagonalization(A,B,max_iter,p,reorth);
 if breaks, max_iter = breaks; end
 [W_T,LambdaT]=eig(T_kk(1:p*max_iter,1:p*max_iter));
 lambdaT=diag(LambdaT);
 [sa vs] = sort(abs(lambdaT),'descend');
 W_T = W_T(:,vs);
 lambdaT= lambdaT(vs);
 R_lambda = zeros(max_iter,1);
 vk = [1:max_iter]';
 for k = 1:max_iter
     k1=round(p*k/3);
     if k1==0
         k1=1;
     end
     k1;
     xx=zeros(k1,1);
     for i=1:k1
         xx(i,1)=abs(lambdaT(i,1)-lambdaA(i,1))/abs(lambdaA(i,1));
     end
    R_lambda(k,1)=max(xx);
 end
  figure(1)
 semilogy(vk, R_lambda, 'o-')
 set(gca,'fontsize',12)
 xlim(vk([1 end]))
 %ylim(ylimit)
 title(titolo);
 grid
% 

%%% Aziz ends here 

% [Q_k,Q_kk,T_kk,rhs,C,normB,breaks]= Block_Lanczos_tridiagonalization(A,B,max_iter,p,reorth);
% if breaks, max_iter = breaks; end
% [W_T,LambdaT]=eig(T_kk(1:p*max_iter,1:p*max_iter));
% lambdaT=diag(LambdaT);
% [sa vs] = sort(abs(lambdaT),'descend');
% W_T = W_T(:,vs);
% lambdaT= lambdaT(vs);
% 
% R_lambda = zeros(p*max_iter,1);
% 
% vk = [1:p*max_iter]';
% for k = 1:p*max_iter
%     k1=round(k/3);
%     if k1==0
%         k1=1;
%     end
%     k1;
%     xx=zeros(k1,1);
%     for j=1:k1
%         xx(j,1)=abs(lambdaT(j,1)-lambdaA(j,1))/abs(lambdaA(j,1));
%     end
%    R_lambda(k,1)=max(xx);
% end
% 
% 
% figure(1)
% semilogy(vk, R_lambda, 'o-')
% set(gca,'fontsize',12)
% xlim(vk([1 end]))
% %ylim(ylimit)
% title(titolo);
% grid
