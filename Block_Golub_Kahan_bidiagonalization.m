% This function to apply m steps of the block Golub—Kahan bidiagonalization
% Input:
% - A non-symmetric matrix A n by l
% - B rhs of the available sys " block vector" n by p
% - p size of the block
% - m number of iterations we want to apply
% Output
% - Q_m l by pm orthogonal
% - Q_mm l by p(m+1) orthogonal
% - P_m n by pm orthogonal
% - P_mm n by p(m+1) orthogonal
% - C_m pm by pm lower block bidiagonal;
% - C_mm p(m+1) by pm lower  block bidiagonal;
% - the rhs of the reduced problem
function [Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks]= Block_Golub_Kahan_bidiagonalization(A,B,m,p,reorth)

if nargin<5, reorth = 1; end

[n,l]=size(A);
sqrteps=sqrt(eps);
Q=zeros(l,p,m);   % initialize Q
P=zeros(n,p,m+1);   % initialize P
L=zeros(p,p,m);   % initialize L
R=zeros(p,p,m+1);   % initialize R
[P(:,:,1),R(:,:,1)]=qr(B,0);
E_1=zeros(p*(m+1),p);
for i=1:p
	E_1(i,i)=1;
end
rhs=E_1*R(:,:,1);
Z=A'*P(:,:,1);
[Q(:,:,1),S]=qr(Z,0);
L(:,:,1)=S';
Q(:,:,1)=Q(:,:,1)/norm(Q(:,:,1));
breaks = 0;
j = 0;
while (j<m) && ~breaks
	j = j+1;
	W=A*Q(:,:,j)-P(:,:,j)*L(:,:,j);
	check=norm(W);
	tau=1e-12;
	if (j*p<n) && (check < tau)
		warning(sprintf('breakdown happened at it. %d',j));
		breaks = j;
	end
	[P(:,:,j+1),R(:,:,j+1)]=qr(W,0);
	if reorth	% reorthognalization P(:,:,j+1)
		for i=1:j
			C = P(:,:,i)'*P(:,:,j+1);
			P(:,:,j+1)=P(:,:,j+1)-P(:,:,i)*C;
			nr = 0;
			while (norm(C) > sqrteps) && (nr<10)
				nr = nr+1;
				C = P(:,:,i)'*P(:,:,j+1);
				P(:,:,j+1) = P(:,:,j+1)-P(:,:,i)*C;
			end
			if nr > 1	
				warning(sprintf('%d reorth steps at it. %d',nr+1,j))
			end
		end
		% normalize
		P(:,:,j+1)=P(:,:,j+1)/norm(P(:,:,j+1));
	end
	if j<m
		Z=A'*P(:,:,j+1)-Q(:,:,j)*R(:,:,j+1)';
		check1=norm(Z);
		tau=1e-12;
		if check1 < tau
			warning(sprintf('breakdown happened at it. %d',j));
			breaks = j;
		end
		[Q(:,:,j+1),S]=qr(Z,0);
		L(:,:,j+1)=S';
		if reorth	% reorthognalization Q(:,:,j+1)
			for i=1:j
				C = Q(:,:,i)'*Q(:,:,j+1);
				Q(:,:,j+1)=Q(:,:,j+1)-Q(:,:,i)*C;
				nr = 0;
				while (norm(C) > sqrteps) && (nr<10)
					nr = nr+1;
					C = Q(:,:,i)'*Q(:,:,j+1);
					Q(:,:,j+1) = Q(:,:,j+1)-Q(:,:,i)*C;
				end
				if nr > 1	
					warning(sprintf('%d reorth steps at it. %d',nr+1,j))
				end
			end
			% normalize
			Q(:,:,j+1)=Q(:,:,j+1)/norm(Q(:,:,j+1));
		end
	end
end
if breaks, m = breaks; end
%%% construct the matrices
for i=1:m
	P_m(:,((i-1)*p)+1:i*p)=P(:,:,i);
end
for i=1:m+1
	P_mm(:,((i-1)*p)+1:i*p)=P(:,:,i);
	%Q_mm(:,((i-1)*p)+1:i*p)=Q(:,:,i);
end
for i=1:m
	Q_m(:,((i-1)*p)+1:i*p)=Q(:,:,i);
	D(((i-1)*p)+1:i*p,((i-1)*p)+1:i*p)=L(:,:,i);
end
LSupd=zeros(p*m);
for i=2:m+1
	LSupd(((i-1)*p)+1:i*p,((i-2)*p)+1:(i-1)*p)=R(:,:,i);
end
C_mm=zeros(p*(m+1),p*m);
C_mm(p*m+1:p*m+p,p*m-(p-1):p*m)=R(:,:,m+1);
C_mm(1:p*m,1:1:p*m)=D+LSupd(1:p*m,1:p*m);
C_m=C_mm(1:p*m,1:p*m);
%%
[s]=svd(A);
RHS=zeros(m,1);
LHS=zeros(m,1);
%    product=L(:,:,2)'*R(:,:,2);
% LHS(1,1)=norm(product);
% RHS(1,1)=s(1,1).^2;
%     for i=2:m-1
%         i
%        RHS(i,1)=s(i,1).^2*RHS(i-1,1);
%
%        product=L(:,:,i+1)'*R(:,:,i+1)*product;
%             LHS(i,1)=norm(product);
%
%     end
%     LHS
%     vk = [1:m]';
%     RHS=RHS./RHS(1,1);
%    figure(1); semilogy(vk,RHS,'o-b',vk, LHS,'+--r');
%   grid on;
% title('Heat');
% TEST=[LHS,RHS]
end
%
