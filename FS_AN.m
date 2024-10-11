function [W,obj,S,p, fn]=FS_AN(X,c,feature_num,k1,beta,k)
%X=X^(dxn)

fea = X';
[nSmp,nFea] = size(fea);
%  options = [];
%  options.NeighborMode = 'KNN';
%  options.k = k;
%   options.NeighborMode = 'Supervised';
%   options.gnd = gnd;
%  options.WeightMode = 'HeatKernel';
%  options.t = 1;
% S = constructW(fea,options); 
% k = 15;
num = nSmp;
distX = L2_distance_1(X,X);
%distX = sqrt(distX);
[distX1, idx] = sort(distX,2);
S = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
% if gamma <= 0
%     gamma = mean(rr);
% end;
lambda = mean(rr);
% T = S;
% P = eye(nSmp);
P=rand(nSmp,1);
P=diag(P./sum(P)); 
% q = P;
% Y1 = zeros(nSmp,1);
% Y2 = zeros(nSmp);
% y3 = 0;
% y4 = 0;
% mu = 0.1;
% r = 1.01;
S0 = (S+S')/2;
D0 = diag(sum(S0));
L0 = D0 - S0;
[F, temp, evs]=eig1(L0, c, 0);

H = eye(num)-1/num*ones(num);
St = X*H*X';
    [U1,S1,V1]=svd(St,'econ');
    T1 = S1;
    T1(S1~=0)=1./S1(S1~=0);
    InvSt = V1*T1'*U1';
% invSt = inv(St);
M = (X*L0*X');
W = eig1(M, feature_num, 0, 0);
% St = X*X';

iter=1;
err=1;
fn = zeros(50, 1);
obj = zeros(50,1);
while err>1e-3 && iter <= 50
        
     distf = L2_distance_1(F',F');
    distx = L2_distance_1(W'*X,W'*X);
%     if iter>5
%         [temp, idx] = sort(distx,2);
%     end;
% %     S = zeros(nSmp);
%     for i=1:nSmp
% %         if islocal == 1
% %             idxa0 = idx(i,2:k+1);
% %         else
%             idxa0 = 1:nSmp;
% %         end;
%         dfi = distf(i,idxa0);
%         dxi = distx(i,idxa0);
%         ad = -(dxi*P+beta*dxi+lambda*dfi)/(2*gamma*beta);
%         S(i,idxa0) = EProjSimplex_new(ad);
%     end;
%     S = S*P;
g = diag(distx*S);
    [p,alpha]=RSWL(g,k1);

    [S,phi]=RSWL(distx*P+distx.*beta+lambda*distf,k);
    P = diag(p);
    S  = (S+S')./2;
    D = diag(sum(S));
    L = D-S;
    
    A = (diag(sum(P*S))+diag(sum(P*S,2))-2*P*S)...
        +beta.*(diag(sum(S))+diag(sum(S,2))-2*S);
    M = InvSt*(X*A*X');
    M = real(M);
    M(isnan(M)==1)=0;
    W = eig1(M,feature_num,0,0);
%     g = -(sum(distx*S,2)./(2*alpha));
%     p = EProjSimplex_new(g);
        
        obj(iter)=sum(P*g)+alpha.*sum(p.^2)...
        +beta*(sum(sum(distx.*S))+sum(phi.*sum(S.^2,2)))...
        +sum(sum(2*lambda*distf.*S));
      
    if iter>1
    err=abs(obj(iter)-obj(iter-1));
    end
    
     F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;

    %eigGap(iter) = ev(c+1)- ev(c);

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    fn(iter) = ev(c+1) - ev(c);
    if fn1 > 0.000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    elseif iter >= 5
         break;
    end

iter=iter+1;    
end
disp(['convergence is',num2str(iter)]);
    
    
    
    
    

function [w,lam]=RSWL(f,k) %k<N
[N,d]=size(f);
w=zeros(N,d);
if d>1
[P,~]=sort(f,2,'ascend');
w=max(((P(:,k+1)*(ones(1,N))-f))./(((k*P(:,k+1)-sum(P(:,1:k),2)))*(ones(1,N))),0);
lam=(k*P(:,k+1)-sum(P(:,1:k),2))/2;
% w = w - diag(diag(w));
else
 [P,~]=sort(f,'ascend');
 w=max(((P(k+1)*ones(N,1)-f)/(k*P(k+1)-sum(P(1:k)))),0);
lam=(k*P(k+1)-sum(P(1:k)))/2;
end
