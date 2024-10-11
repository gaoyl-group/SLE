clc;clear;close all;

load("breast.mat")

%% 调试参数
num = size(X,1);
c = length(unique(Y));
dim = c-1;
ps = [1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3];

acc_m = zeros(5,7,5);
acc_r = zeros(5,5);

indices = crossvalind('Kfold', num, 5);
for j = 1:5
    test = (indices == j); train = ~test;
    train_index = find(train==1)';
    test_index = find(test==1)';
    X_train = X(train_index,:);
    train_label = Y(train_index);
    X_test = X(test_index,:);
    test_label = Y(test_index);

    num2 = size(X_train,1);
    ps2 = [floor(num2*0.5) floor(num2*0.6) floor(num2*0.7) floor(num2*0.8) floor(num2*0.9)];
    
    %% RPCA-AN
    % australian k = 0.8*num
    % vote k = 0.8*num
    % monk1 k = 0.9*num
    % breast k = 0.7*num
    for ii = 1:5
        [W, ~] = RPCA_AN(X_train', ps2(ii), dim);
        X1 = real(X_train*W);
        X2 = real(X_test*W);
        mdl = ClassificationKNN.fit(X1, train_label,'NumNeighbors',5);
        predict_label = predict(mdl,X2);
        acc_r(ii,j)=length(find(predict_label == test_label))/length(test_label)*100;
    end

    %% OUR
    % australian k = 0.8*num, beta = 1e-3;
    % vote k = 0.9*num, beta = 1e-1
    % monk1 k = 0.8*num beta - 1e2
    % breast k = 0.5*num beta = 1e-2
    for ii = 1:5
        for k = 1:7
            [W,~] = FS_AN(X_train',c,dim,ps2(ii),ps(k),20);
            X1 = real(X_train*W);
            X2 = real(X_test*W);
            mdl = ClassificationKNN.fit(X1, train_label,'NumNeighbors',5);
            predict_label = predict(mdl,X2);
            acc_m(ii,k,j) = length(find(predict_label == test_label))/length(test_label)*100;
        end
    end
end

save('breast_tune_p.mat', 'acc_m', 'acc_r');

