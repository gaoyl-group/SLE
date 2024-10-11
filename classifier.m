clc;clear;close all;

load("breast.mat")
num = size(X,1);
c = length(unique(Y));
dim = c-1;

%% pca
H = eye(num)-1/num*ones(num);
St = X'*H;
[U, S, V] = svd(St,'econ'); s = diag(S);
Wp = U(:,1:dim);
X_pca = X*Wp;

%% LPP
H = eye(num)-1/num*ones(num);
A = selftuning(X, 20);
L = diag(sum(A,2))-A;
St =X'*(diag(sum(A,2)))*X;
invSt = pinv(St);
Sl = X'*L*X;
D = diag(sum(A,2));
M = (invSt)*Sl;
[Wl, temp, ev]=eig1(M, dim, 0, 0);
Wl = Wl*diag(1./sqrt(diag(Wl'*Wl)));
X_lpp = real(X*Wl);

%% PCAN
[W, ~] = PCANK(X',c,dim,20,-1,0);
X_PCAN = real(X*W);

%% RPCA-AN
[W, ~] = RPCA_AN(X', floor(0.7*num), dim);
X_RPCA = real(X*W);

%% OUR
[W,~] = FS_AN(X',c,dim,floor(0.5*num),1e-2,20);
X_OUR = real(X*W);

save('breast_reduce_data.mat', 'X', 'Y', 'X_pca', 'X_lpp', 'X_PCAN', ...
    'X_RPCA', 'X_OUR')