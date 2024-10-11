function [W,WResult] = InterationW_SPCA(X,gamma,m)
%X: data matrix(dim*num)
%gamma: regularization parameter
%m: projection dimension of W (dim*m)

num = size(X,2);
dim = size(X,1);

INTER_W = 30;
Q = eye(dim);

H = eye(num)-(1/num)*ones(num);
St = X*H*X';
St = -St; 

p=1; % L_2p
for i = 1:INTER_W
    
    tempStQ = (St+gamma*Q);
    [vec,val] = eig(tempStQ);
    [~,di] = sort(diag(val));
    W = vec(:,di(1:m));

    tempQ = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    Q = diag(tempQ);

    w1(i) = trace(W'*St*W); %  Tr(W'*St*W)
    w2(i) = gamma*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
    WResult(i) = w1(i)+w2(i);
    disp(['iter', num2str(i)]);
    if i > 1 && abs(WResult(i-1)-WResult(i)) < 1e-3
        break;
    end;
    
end;
end