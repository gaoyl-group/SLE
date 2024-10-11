function [W, A, F, G, y, ObjHistory] = DCLA(X, c, d, k, r, lambda, islocal)
% X: dim*num data matrix, each column is a data point
% c: number of clusters
% d: projected dimension
% k: number of neighbors to determine the initial graph, and the parameter r if r<=0
% r: paremeter,. If r<0, then it is determined by algorithm with k
% islocal: 
%           1: only update the similarities of the k neighbor pairs, the neighbor pairs are determined by the distances in the original space 
%           0: update all the similarities
% W: dim*d projection matrix
% y: num*1 cluster indicator vector
% A: num*num learned symmetric similarity matrix
% F: num*c cluster indicator matrix
% G: d*c cluster centroid matrix

NITER = 30;
[dim, num] = size(X);
if nargin < 7
    islocal = 1;
end
if nargin < 6
    lambda = 1e-4;
end
if nargin < 5
    r = -1;
end
if nargin < 4
    k = 10;
end
if nargin < 3
    d = c-1;
end


%% initial similarity matrix, r, lambda, projected matrix
% initial W
%W = orth(rand(dim,d));
coeff = pca(X'); 
W = orth(coeff);
W = W(:,1:d);
% initial A,r
distX = L2_distance_1(W'*X, W'*X);
[distX1, idx] = sort(distX, 2);
A = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i, 2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i, 2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
if r < 0
    r = mean(rr);
end

A = (A+A')/2;
D = diag(sum(A));
L = D - A;

% initial F
F = zeros(num, c);
col_F = randi(c, 1, num);
ind_F = sub2ind([num,c], 1:num, col_F);
F(ind_F) = 1;

% initial G
temp_dim = randperm(dim);
temp_num = randperm(num);
G = X(temp_dim(1:d), temp_num(1:c));
% calcultater obj
ObjHistory = zeros(NITER,1);
%obj = sum(sum(distwx.*S0))+gamma*sum(sum(S0.^2))+lamda*norm( W'*X-G*F' ,'fro');

%% iterative algorithm
for iter = 1:NITER
    %ObjHistory(iter) = obj;
    % update F
    F_better = F;
    objKM_better = norm((X'*W-F*G'), 'fro')^2;
    find_better = 0;
    for i = 1:10
        col_F = randi(c, 1, num);
        ind_F = sub2ind([num,c], 1:num, col_F);
        F(ind_F) = 1;
        objKM = norm((X'*W-F*G'), 'fro')^2;
        if objKM < objKM_better
            F_better = F;
            find_better = 1;
        end
    end
    F = F_better;
    %if not
    if find_better == 0
        for i = 1:num
            id_F = 1;
            objKM_better = norm((W'*X(:,i)-G(:,1)),'fro')^2;
            for j = 2:c
                objKM = norm((W'*X(:,i)-G(:,j)),'fro')^2;
                if objKM < objKM_better
                    id_F = j;
                    objKM_better = objKM;
                end
            end
            F(i,:) = 0;
            F(i, id_F) = 1;
        end
    end

    %update G
    invF = pinv(F'*F+eps);
    G = W'*X*F*invF;

    %update W
    M1 = L+lambda*eye(num)-lambda*(F*invF*F');
    M = X*M1*X';
    M = (M+M')/2;
    W = eig1(M, d, 0, 0);
    W = W*diag(1./sqrt(diag(W'*W)));

    %update A
    distx = L2_distance_1(W'*X, W'*X);
    if iter>5
        [~,idx] = sort(distx, 2);
    end
    A = zeros(num);
    for i = 1:num
        if islocal == 1
            idxa0 = idx(i, 2:k+1);
        else
            idxa0 = 1:num;
        end
        dxi = distx(i, idxa0);
        ad = -dxi/(2*r);
        A(i,idxa0) = EProjSimplex_new(ad);
    end       

    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;

    %obj = trace(W'*X*L*X'*W)+r*norm(A,'fro')^2+lambda*norm((X'*W-F*G'), 'fro')^2;
    ObjHistory(iter) = sum(sum(distx.*A))+r*sum(sum(A.^2))+lambda*norm( W'*X-G*F' ,'fro');
    disp(['iter= ' num2str(iter), ',FunVal=' num2str(ObjHistory(iter)) ]);
    if iter>1 && abs(ObjHistory(iter)-ObjHistory(iter-1)) < 1e-10 && ObjHistory(iter) < ObjHistory(iter-1)
        break;
    end
end
disp(['iter= ' num2str(iter), ',FunVal=' num2str(ObjHistory(iter)) ]);
%% calculate label
%y = vec2ind(F')';
y = make_label(F);
end

