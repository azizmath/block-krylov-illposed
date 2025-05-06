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
function [Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks]= Block_Golub_Kahan_bidiagonalization_IRtools(A,B,m,p)
    n = size(B,1);
    P=zeros(n,p,m+1);   % initialize P 
    R=zeros(p,p,m+1);   % initialize R 
    [P(:,:,1),R(:,:,1)]=qr(B,0);
    Z = Atransp_times_vec(A,P(:,:,1));
    l = size(Z,1);
    sqrteps=sqrt(eps);
    Q=zeros(l,p,m);   % initialize Q 
    L=zeros(p,p,m);   % initialize L 
    
    E_1=zeros(p*(m+1),p);
    
    for i=1:p
        E_1(i,i)=1;
    end
    rhs=E_1*R(:,:,1);
    
    [Q(:,:,1),S]=qr(Z,0);
    L(:,:,1)=S';
    Q(:,:,1)=Q(:,:,1)/norm(Q(:,:,1));
    for j=1:m
        W=A_times_vec(A,Q(:,:,j))-P(:,:,j)*L(:,:,j);
        check=norm(W);
        tau=1e-12;
        if check < tau
            warning('breaks down may happen',j);
            breaks=0;
        else 
            breaks=1;
        end
        [P(:,:,j+1),R(:,:,j+1)]=qr(W,0);
        % Full reorthogonalization step.
        for i=1:j
            C=(P(:,:,j+1)'*P(:,:,i))';
            P(:,:,j+1)=P(:,:,j+1)-P(:,:,i)*C;
            if norm(C)>sqrteps
                C2 =(P(:,:,j+1)'*P(:,:,i))'; 
                P(:,:,j+1) = P(:,:,j+1)-P(:,:,i)*C2; 
            end
        end
        if j<m
            Z=Atransp_times_vec(A,P(:,:,j+1))-Q(:,:,j)*R(:,:,j+1)';
            check1=norm(Z);
             tau=1e-12;
            if check1 < tau
                warning('breaks down may happen',j);
                breaks=0;
            else 
                breaks=1;
            end
            [Q(:,:,j+1),S]=qr(Z,0);
            L(:,:,j+1)=S';
            % reorthognalization Q(:,:,j+1)
             for i=1:j
                C3=(Q(:,:,j+1)'*Q(:,:,i))'; 
                Q(:,:,j+1)=Q(:,:,j+1)-Q(:,:,i)*C3;      
                if norm(C3)>sqrteps
                    C4 =(Q(:,:,j+1)'*Q(:,:,i))'; 
                     Q(:,:,j+1) = Q(:,:,j+1)-Q(:,:,i)*C4; 
                end
             end    
            % normalize
            Q(:,:,j+1)=Q(:,:,j+1)/norm(Q(:,:,j+1));
        end
    end
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
%     [s]=svd(A);
%     RHS=zeros(m,1);
%     LHS=zeros(m,1);
%     product=L(:,:,2)'*R(:,:,2);
%     LHS(1,1)=norm(product);
%     RHS(1,1)=s(1,1).^2;;
%     for i=2:m-1
%         i
%        RHS(i,1)=s(i,1).^2*RHS(i-1,1);
%        product=L(:,:,i+1)'*R(:,:,i+1)*product;
%             LHS(i,1)=norm(product);
%         
%     end
%     LHS
%     vk = [1:m]';
%    figure(1); semilogy(vk,RHS,'o-b',vk, LHS,'+--r'); 
%   grid on;
% title('Heat'); 
end
% 
