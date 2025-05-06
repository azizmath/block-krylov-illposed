% This code provides us with the results for  plots of the figures
% (Fig4.4, Fig4.5, Fig4.6, and Fig4.7)
clear
n=200;p=5;noiselevel=10^-2;max_iter=40;
ylimit = [1e-16 1e1];
makefigs = 0;
reorth = 1;

%[A] = deriv2(n,1); [AA,b,x] = phillips(n);A=-A; titolo = 'Deriv2';
%[A,b,x] = foxgood(n); titolo = 'Foxgood';
%[A] = gravity(n,1,0,1,0.25); [AA,b,x] = phillips(n); titolo = 'Gravity';
[A,b,x] = phillips(n); titolo = 'Phillips';

%[A,b,x] = shaw(n); titolo = 'Shaw';
%[A,b,x] = heat(n,1); titolo = 'Heat';
%A = gallery('lotkin',n);[At b x] = shaw(n); titolo = 'Lotkin';

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
R_wk=zeros(max_iter,1);

[Q_k,Q_kk,T_kk,rhs,C,normB,breaks]= Block_Lanczos_tridiagonalization(A,B,max_iter,p,reorth);
if breaks, max_iter = breaks; end

% first test
if breaks
	k = breaks;
else
	k = round(max_iter/2);
end
vi = [1:p*k]';
ind1 = zeros(p*k,1);
P = eye(n) - Q_k(:,1:p*k)*Q_k(:,1:p*k)';
for i = 1:p*k
	ind1(i) = norm(P*W(:,1:i),'fro');
end

figure(1)
semilogy(vi, ind1, 'o-')
set(gca,'fontsize',12)
xlim(vi([1 end]))
ylim(ylimit)
title(titolo);
grid
if makefigs
	filnam = [titolo 'L1.eps'];
	print('-depsc2',filnam)
end

% second and third test
ind2 = zeros(max_iter,1);
ind3 = zeros(max_iter,1);
vk = [1:max_iter]';
for k = 1:max_iter
	P = eye(n) - Q_k(:,1:p*k)*Q_k(:,1:p*k)';
	k3 = round(p*k/3);
	ind2(k) = norm(P*W(:,1:k3),'fro');
	k2 = round(p*k/2);
	ind3(k) = norm(P*W(:,1:k2),'fro');
end

figure(2)
semilogy(vk, ind2, 'o-')
set(gca,'fontsize',12)
xlim(vk([1 end]))
ylim(ylimit)
title(titolo);
grid
if makefigs
	filnam = [titolo 'L2.eps'];
	print('-depsc2',filnam)
end

figure(3)
semilogy(vk, ind3, 'o-')
set(gca,'fontsize',12)
xlim(vk([1 end]))
ylim(ylimit)
title(titolo);
grid
if makefigs
	filnam = [titolo 'L3.eps'];
	print('-depsc2',filnam)
end

