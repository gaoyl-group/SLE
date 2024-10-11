clc;clear;close all;

addpath("gefuns\")

[X,ColorVector] = Generate3DClustersData();
cc = ColorVector;
X = X';
[ M , N ] = size(X);
X_mean = mean( X , 2 );
X = X - repmat( X_mean , 1 , N ); 
X = X';
num = N;
c = 3;
%% LPP
H = eye(num)-1/num*ones(num);
A = selftuning(X, 20);
L = diag(sum(A,2))-A;
St =X'*(diag(sum(A,2)))*X;
invSt = pinv(St);
Sl = X'*L*X;
D = diag(sum(A,2));
M = (invSt)*Sl;
[Wl, temp, ev]=eig1(M, 2, 0, 0);
Wl = Wl*diag(1./sqrt(diag(Wl'*Wl)));
Y1 = real(X*Wl);

%% PCAN
[W, ~] = PCANK(X',c,2,20,-1,0);
Y2 = real(X*W);

%% OUR
[W,~] = FS_AN(X',c,2,floor(0.85*N),0.01,20);
Y3 = real(X*W);

%% PLOT
dot = 10;
C = X;
set(figure,'position',[200,200,300,300] );set(gca,'looseInset',[0 0 0 0]); colormap(jet);
scatter3( C(:,1) , C(:,2) , C(:,3) , dot , cc , 'filled' );title('original')
set(figure,'position',[200,200,300,300] );set(gca,'looseInset',[0 0 0 0]); colormap(jet);
scatter( Y1(:,1) , Y1(:,2) , dot , cc , 'filled');title('LPP')
set(figure,'position',[200,200,300,300] );set(gca,'looseInset',[0 0 0 0]); colormap(jet);
scatter( Y2(:,1) , Y2(:,2) , dot , cc , 'filled');title('PCAN')
set(figure,'position',[200,200,300,300] );set(gca,'looseInset',[0 0 0 0]); colormap(jet);
scatter( Y3(:,1) , Y3(:,2) , dot , cc , 'filled');title('Ours')