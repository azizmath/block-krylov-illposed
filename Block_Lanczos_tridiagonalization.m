% This function to compute m steps of the block lanczos tridiagonalization
% Input:
% - A: symmetric matrix A n by n
% - B: rhs of the availible sys " block vector" n by p
% - p: size of the block
% - m: number of iterations we want to apply
% Output
% - Q_m: n by pm orhtogonal
% - Q_m+1: n by p(m+1) orthognal
% - T_m+1,m: p(m+1) by pm block tridiagonal;
% - rhs: the rhs of the reduced problem
% - C: subdigonal blocks B_i i:1,.....,m+1
% - normB= % vector with norm of the subdiagonal blocks \|B_i\| i:1,.....,m+1
function [Q_m,Q_mm,T_mm,rhs,C,normB,breaks]= Block_Lanczos_tridiagonalization(A,B,m,p,reorth)

if nargin<5, reorth = 1; end

n = size(A,1);
sqrteps=sqrt(eps);
C=zeros(p,p,m+1);
normB=zeros(m,1);
breaks=0;
k = 0;
while (k<m) && ~breaks
	k = k+1;
	[X(:,:,1),S]=qr(B,0);
	S=S;
	E_1=zeros(p*(k+1),p);
	for i=1:p
		E_1(i,i)=1;
	end
	rhs=E_1*S;
	M(:,:,1)=X(:,:,1)'*A*X(:,:,1);
	if k==1
		R=A*X(:,:,1)-X(:,:,1)*M(:,:,1);
		[X(:,:,2),C(:,:,1)]=qr(R,0);
		normB(1,1)=norm(C(:,:,1));
		X(:,:,2)=X(:,:,2)-X(:,:,1)*(X(:,:,1)'*X(:,:,2));
		M(:,:,2)=X(:,:,2)'*A*X(:,:,2);
	else
		R=A*X(:,:,k)-X(:,:,k)*M(:,:,k)-X(:,:,k-1)*C(:,:,k-1)';
		check=norm(R); % check breakdown
		tau=1e-12;
		if (k*p<n) && (check < tau)
			warning(sprintf('breakdown happened at it. %d',k));
			breaks=k;
		else
			breaks=0;
		end
	end
	[X(:,:,k+1),C(:,:,k)]=qr(R,0);
	normB(k,1)=norm(C(:,:,k));
	if reorth
		for j=1:k % reorthogonalization
			H = X(:,:,j)'*X(:,:,k+1);
			X(:,:,k+1)=X(:,:,k+1)-X(:,:,j)*H;
			nr = 0;
			while (norm(H) > sqrteps) && (nr<10)
				nr = nr+1;
				H = X(:,:,j)'*X(:,:,k+1);
				X(:,:,k+1) = X(:,:,k+1)-X(:,:,j)*H;
			end
			if nr > 1	
				warning(sprintf('%d reorth steps at it. %d',nr+1,j))
			end
		end
		X(:,:,k+1)=X(:,:,k+1)/norm(X(:,:,k+1));
	end
	M(:,:,k+1)=X(:,:,k+1)'*A*X(:,:,k+1);
	D=zeros(p*(k+1),p*k); %  consruct matrix Q_m and matrix T_m1
	for i=1:k
		Q_m(:,((i-1)*p)+1:i*p)=X(:,:,i);
		D(((i-1)*p)+1:i*p,((i-1)*p)+1:i*p)=M(:,:,i);
	end
	for i=1:k+1
		Q_mm(:,((i-1)*p)+1:i*p)=X(:,:,i);
	end
	LSupd=zeros(p*(k+1),p*k);
	USupd=zeros(p*(k+1),p*k);
	for i=2:k+1
		LSupd(((i-1)*p)+1:i*p,((i-2)*p)+1:(i-1)*p)=C(:,:,i-1);
	end
	for i=2:k
		USupd(((i-2)*p)+1:(i-1)*p,((i-1)*p)+1:i*p)=C(:,:,i-1)';
	end
	T_mm=D+LSupd+USupd;
end
%     fff=zeros(m,1);
% fff(1,1)=norm(A*Q_m(:,1:p));
% for i=2:m
%     fff(i,1)=norm(A*Q_mm(:,p+(i-1):i*p));
% end
% figure (1); semilogy(fff,'b-o'); title ('Phillips');
end
