% This code provides us with the results for  plots of the figures
% (Fig4.4, Fig4.5, Fig4.6, and Fig4.7)
clear
n=400;p=5;noiselevel=0;max_iter=n/p;
ylimit = [1e-16 1e1];
makefigs = 0;
reorth = 1;

%[A,b,x] = baart(n); titolo = 'Baart';
%[A,b,x] = heat(n,1); titolo = 'Heat';
%A = gallery('lotkin',n);[At, b, x] = shaw(n); titolo = 'Lotkin';
%[A,b,x] = wing(n,1/3,2/3); titolo = 'Wing';
%N=20;[A,b,x] = tomo(N,1);n=N^2;titolo = 'tomo';
[A,b,x] = fulltomo(n);titolo = 'tomo';
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
% 
[W,LambdaA]=eig(A);
lambdaA=diag(LambdaA);
[sa vs] = sort(abs(lambdaA),'descend');
W = W(:,vs);
lambdaA= lambdaA(vs);

[U,Sigma,V]=csvd(A);

[Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks]= Block_Golub_Kahan_bidiagonalization(A,B,max_iter,p,reorth);

if breaks, max_iter = breaks; end

% first test
if breaks
	k = breaks;
else
	k = round(max_iter/2);
end
vi = [1:p*k]';
ind1 = zeros(k,1);
P = eye(n) - Q_m(:,1:p*k)*Q_m(:,1:p*k)';
Q = eye(n) - P_m(:,1:p*k)*P_m(:,1:p*k)';
for i = 1:p*k
	aa = norm(P*V(:,1:i),'fro');
	bb = norm(Q*U(:,1:i),'fro');
	ind1(i) = max(aa,bb);
end

figure(1)
semilogy(vi, ind1, 'o-')
set(gca,'fontsize',12)
xlim(vi([1 end]))
ylim(ylimit)
title(titolo);
grid
if makefigs
	filnam = [titolo 'GK1.eps'];
	print('-depsc2',filnam)
end

% second and third test
ind2 = zeros(max_iter,1);
ind3 = zeros(max_iter,1);
vk = [1:max_iter]';
for k = 1:max_iter
	P = eye(n) - Q_m(:,1:p*k)*Q_m(:,1:p*k)';
	Q = eye(n) - P_m(:,1:p*k)*P_m(:,1:p*k)';
	k3 = round(p*k/3);
	aa = norm(P*V(:,1:k3),'fro');
	bb = norm(Q*U(:,1:k3),'fro');
	ind2(k) = max(aa,bb);
	k2 = round(p*k/2);
	aa = norm(P*V(:,1:k2),'fro');
	bb = norm(Q*U(:,1:k2),'fro');
	ind3(k) = max(aa,bb);
end

figure(2)
semilogy(vk, ind2, 'o-')
set(gca,'fontsize',12)
xlim(vk([1 end]))
ylim(ylimit)
title(titolo);
grid
if makefigs
	filnam = [titolo 'GK2.eps'];
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
	filnam = [titolo 'GK3.eps'];
	print('-depsc2',filnam)
end

