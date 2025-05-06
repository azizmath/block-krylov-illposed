%clear ;clf;close all;clc;
n = 256;
opt = PRset('PSF', 'shake'); A = PRblur(n, opt);
rng(0)

Kall = kronApprox(A);
A1 = Kall.a{1};
B1 = Kall.b{1};

X = imread('tissue.png');
X = X(125+(1:n),1:n,:);
X = im2double(X);

%figure, imshow(X, [])

XX = [reshape(X(:,:,1), n^2,1), reshape(X(:,:,2), n^2,1), reshape(X(:,:,3), n^2,1)];
A = @(xx,tflag) OPkroncol(xx,B1,A1,tflag);
BB = A_times_vec(A,XX);
[bn, NoiseInfo] = PRnoise(BB);
nrmbn = norm(bn, 'fro');

%figure, imshow(reshape(bn, n, n, 3), [])

% solution via block method
max_iter = 50;
p = 3;
[~,R(:,:,1)] = qr(bn,0);
[Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks] = Block_Golub_Kahan_bidiagonalization_IRtools(A,bn,max_iter,p);
Enrm_BGKB = zeros(max_iter,1);
Rnrm_BGKB = zeros(max_iter,1);
for i=1:max_iter
    E_1=zeros(p*(i+1),p);
    E_1(1:p,1:p) = eye(p);
    rhs=E_1*R(:,:,1);
    [U,Sigma,V]=svd(C_mm(1:(i+1)*p,1:i*p),'econ'); % Tik
    ss=diag(Sigma);
    mu = 0; % no regularization
    beta2 = ss./(ss.^2+mu.^2);
    Y_mu= V*diag(beta2)*U'*rhs; X_m=Q_m(:,1:i*p)*Y_mu;% Tik
    Enrm_BGKB(i) = norm(X_m-XX,'fro')/norm(XX,'fro');
    Rnrm_BGKB(i) = norm(rhs - C_mm(1:(i+1)*p,1:i*p)*Y_mu,'fro')/norm(bn,'fro');
end

% solution via block method + Tikhonov regularization
norm_error = 1e-2*nrmbn;
tau = 1.1;
Enrm_BGKBtik = zeros(max_iter,1);
Rnrm_BGKBtik = zeros(max_iter,1);
RegP = zeros(max_iter,1);
[~,R(:,:,1)] = qr(bn,0);
for tries = 1:8
tic
[Q_m,P_m,P_mm,C_m,C_mm,rhs,breaks] = Block_Golub_Kahan_bidiagonalization_IRtools(A,bn,max_iter,p);
for i=1:max_iter
    E_1=zeros(p*(i+1),p);
    E_1(1:p,1:p) = eye(p);
    rhs=E_1*R(:,:,1);
    % solve the reduced problem by Tikhonov regularization
    [U,Sigma,V]=svd(C_mm(1:(i+1)*p,1:i*p),'econ'); % Tik
    ss=diag(Sigma);
    mu = 0;
    beta2 = ss./(ss.^2+mu.^2);
    Y_mu= V*diag(beta2)*U'*rhs;
    Rnrm_BGKBtik(i) = norm(rhs - C_mm(1:(i+1)*p,1:i*p)*Y_mu,'fro')/nrmbn;
    if Rnrm_BGKBtik(i)<= tau*1e-2
        mu = fzero(@(l)blockdiscrfcn(l, C_mm(1:(i+1)*p,1:i*p), eye(i*p), rhs, nrmbn, tau*1e-2), [0, 1e10]);
        Rnrm_BGKBtik(i) = norm(rhs - C_mm(1:(i+1)*p,1:i*p)*Y_mu,'fro')/norm(bn,'fro');
    end
    RegP(i) = mu;
    beta2 = ss./(ss.^2+mu.^2);
    Y_mu= V*diag(beta2)*U'*rhs; X_m=Q_m(:,1:i*p)*Y_mu;% Tik
    Enrm_BGKBtik(i) = norm(X_m-XX,'fro')/norm(XX,'fro'); 
end
time_block(tries) = toc;
end

% solution via std lsqr
Avec = @(xx,tflag) OPkroncolvec(xx,B1,A1,tflag);
bnvec = bn(:);
K = 150;
opt = IRset('NoStop', 'on', 'RegParam', 0, 'x_true', XX(:));
[X_lsqr, info_lsqr] = IRhybrid_lsqr(Avec, bnvec, K, opt);

% solution via std lsqr + Tikhonov regularization
Avec = @(xx,tflag) OPkroncolvec(xx,B1,A1,tflag);
bnvec = bn(:);
K = 150;
opt = IRset(opt, 'RegParam', 'discrepit', 'NoiseLevel', 1e-2);
for tries = 1:8
tic
[X_hlsqr, info_hlsqr] = IRhybrid_lsqr(Avec, bnvec, K, opt);
time_std(tries) = toc;
end

% estimate (3.11)
sA1 = svd(A1);
sB1 = svd(B1);
sA = sort(kron(sA1,sB1), 'descend');
%
Prod = eye(p);
lhs_311 = zeros(max_iter-1,1);
for m = 1:max_iter-1
    Prod = C_mm(p*m+1:p*m+p, p*m-(p-1):p*m)*Prod;
    Prod = C_mm(p*m+1:p*m+p,p*m+1:p*m+p)'*Prod;
    lhs_311(m) = norm(Prod);
end
%
rhs_311 = zeros(max_iter-1,1);
for m = 1:max_iter-1
    rhs_311(m) = prod(sA(1:m))^2;
end
