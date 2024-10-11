clc;clear;close all

datasets = ["cars", "vote", "glass", "australian", "lung_discrete", "TOX_171", "Prostate_GE", "GLIOMA", "nci9", "bupa"];

for i = 1:length(datasets)
    data = strcat(datasets(i), '.mat');
    load(data)
    [numX, ~] = size(X);
    c = length(unique(Y));

    nrepeat = 10;
    dim = c - 1;
    indices = crossvalind('Kfold', numX, nrepeat);

    accuracy1 = zeros(1, 10);precision1 = zeros(1, 10);recall1 = zeros(1, 10);F11 = zeros(1, 10);
    accuracy2 = zeros(1, 10);precision2 = zeros(1, 10);recall2 = zeros(1, 10);F12 = zeros(1, 10);
    accuracy3 = zeros(1, 10);precision3 = zeros(1, 10);recall3 = zeros(1, 10);F13 = zeros(1, 10);
    accuracy4 = zeros(1, 10);precision4 = zeros(1, 10);recall4 = zeros(1, 10);F14 = zeros(1, 10);
    accuracy5 = zeros(1, 10);precision5 = zeros(1, 10);recall5 = zeros(1, 10);F15 = zeros(1, 10);
    accuracy6 = zeros(1, 10);precision6 = zeros(1, 10);recall6 = zeros(1, 10);F16 = zeros(1, 10);
    accuracy7 = zeros(1, 10);precision7 = zeros(1, 10);recall7 = zeros(1, 10);F17 = zeros(1, 10);
    for iter = 1:nrepeat
        test = (indices == iter); train = ~test;
        train_index = find(train == 1)';
        test_index = find(test == 1)';
        train_data = X(train_index, :)';train_label = Y(train_index);
        test_data = X(test_index, :)';test_label = Y(test_index);

        %% PCA
        num = size(train_data, 2);
        H = eye(num)-1/num*ones(num);
        St = train_data*H;
        [U, S, V] = svd(St,'econ'); s = diag(S);
        Wp = U(:,1:dim);
        X1 = real(train_data'*Wp); 
        X2 = real(test_data'*Wp);
        mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
        predict_label = predict(mdl, X2);
        [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
        accuracy1(iter) = accuracy;
        precision1(iter) = mean(precision);
        recall1(iter) = mean(recall);
        F11(iter) = mean(F1);
        
        %% LPP
        H = eye(num)-1/num*ones(num);
        A = selftuning(train_data', 10);
        L = diag(sum(A,2)) - A;
        St = train_data*(diag(sum(A,2)))*train_data';
        invSt = pinv(St);
        Sl = train_data*L*train_data';
        D = diag(sum(A,2));
        M = (invSt)*Sl;
        [Wl, temp, ev]=eig1(M, dim, 0, 0);
        Wl = Wl*diag(1./sqrt(diag(Wl'*Wl)));
        X1 = real(train_data'*Wl);
        X2 = real(test_data'*Wl);
        mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
        predict_label = predict(mdl, X2);
        [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
        accuracy2(iter) = accuracy;
        precision2(iter) = mean(precision);
        recall2(iter) = mean(recall);
        F12(iter) = mean(F1);
        
        %% PCAN
        [Wpcan, ~] = PCANK(train_data, c, dim, 10, -1, 0);
        X1 = real(train_data'*Wpcan);
        X2 = real(test_data'*Wpcan);
        mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
        predict_label = predict(mdl, X2);
        [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
        accuracy3(iter) = accuracy;
        precision3(iter) = mean(precision);
        recall3(iter) = mean(recall);
        F13(iter) = mean(F1);

        %% RPCA-AN
        k1 = [floor(num*0.5) floor(num*0.6) floor(num*0.7) floor(num*0.8) floor(num*0.9)];
        acc_p = zeros(1,5);pre_p = zeros(1,5);re_p = zeros(1,5);f_p = zeros(1,5);
        for j = 1:length(k1)
            [Wrpca, ~] = RPCA_AN(train_data, k1(j), dim);
            X1 = real(train_data'*Wrpca);
            X2 = real(test_data'*Wrpca);
            mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
            predict_label = predict(mdl, X2);
            [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
            acc_p(j) = accuracy;
            pre_p(j) = mean(precision);
            re_p(j) = mean(recall);
            f_p(j) = mean(F1);
        end
        accuracy4(iter) = max(acc_p);
        precision4(iter) = max(pre_p);
        recall4(iter) = max(re_p);
        F14(iter) = max(f_p);

        %% DCLA
        [Wdcla, ~] = DCLA(train_data, c, dim, 10, -1, 0);
        X1 = real(train_data'*Wdcla);
        X2 = real(test_data'*Wdcla);
        mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
        predict_label = predict(mdl, X2);
        [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
        accuracy5(iter) = accuracy;
        precision5(iter) = mean(precision);
        recall5(iter) = mean(recall);
        F15(iter) = mean(F1);

        %% SPCAFS
        gamma = [1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3];
        acc_p = zeros(1,7);pre_p = zeros(1,7);re_p = zeros(1,7);f_p = zeros(1,7);
        for j = 1:7
            [id, obj, W] = SPCAFS(train_data, gamma(j), dim);
            X1 = train_data(id(1:dim), :)';
            X2 = test_data(id(1:dim), :)';
            mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
            predict_label = predict(mdl, X2);
            [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
            acc_p(j) = accuracy;
            pre_p(j) = mean(precision);
            re_p(j) = mean(recall);
            f_p(j) = mean(F1);
        end
        accuracy6(iter) = max(acc_p);
        precision6(iter) = max(pre_p);
        recall6(iter) = max(re_p);
        F16(iter) = max(f_p);

        %% SLE
        beta = [1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3];
        k1 = [floor(num*0.5) floor(num*0.6) floor(num*0.7) floor(num*0.8) floor(num*0.9)];
        acc_p = zeros(7,5);pre_p = zeros(7,5);re_p = zeros(7,5);f_p = zeros(7,5);
        for j1 = 1:7
            for j2 = 1:5
                [W, ~] = FS_AN(train_data, c, dim, k1(j2), beta(j1), 10);
                X1 = real(train_data'*W);
                X2 = real(test_data'*W);
                mdl = ClassificationKNN.fit(X1, train_label, 'NumNeighbors', 1);
                predict_label = predict(mdl, X2);
                [confus,accuracy,numcorrect,precision,recall,F1] = compute_accuracy_F (test_label,predict_label,1:c);
                acc_p(j1, j2) = accuracy;
                pre_p(j1, j2) = mean(precision);
                re_p(j1, j2) = mean(recall);
                f_p(j1, j2) = mean(F1);
            end
        end
        accuracy7(iter) = max(acc_p(:));
        precision7(iter) = max(pre_p(:));
        recall7(iter) = max(re_p(:));
        F17(iter) = max(f_p(:));
    end
    file_name = strcat(datasets(i), '_res.mat');
    save(file_name)

    acc = [mean(accuracy1) mean(accuracy2) mean(accuracy3) mean(accuracy4) mean(accuracy5) mean(accuracy6) mean(accuracy7)];
    prec = [mean(precision1) mean(precision2) mean(precision3) mean(precision4) mean(precision5) mean(precision6) mean(precision7)];
    rec = [mean(recall1) mean(recall2) mean(recall3) mean(recall4) mean(recall5) mean(recall6) mean(recall7)];
    f1 = [mean(F11) mean(F12) mean(F13) mean(F14) mean(F15) mean(F16) mean(F17)];

    acc_std = [std(accuracy1) std(accuracy2) std(accuracy3) std(accuracy4) std(accuracy5) std(accuracy6) std(accuracy7)];
    prec_std = [std(precision1) std(precision2) std(precision3) std(precision4) std(precision5) std(precision6) std(precision7)];
    rec_std = [std(recall1) std(recall2) std(recall3) std(recall4) std(recall5) std(recall6) std(recall7)];
    f1_std = [std(F11) std(F12) std(F13) std(F14) std(F15) std(F16) std(F17)];
    disp('[SAVE RESULTS TO TXT]');
    title = ["PCA", "LPP", "PCAN", "RPCAN-AN", "DCLA", "SPCAFS", "OURS"];
    fileName = strcat(datasets(i), '_rest.txt');
    fid = fopen(fileName, 'a');
    fprintf(fid,'%s\r\n',repmat('----', 1, 16)); % 分割线
    fprintf(fid,'%s\t','Measure');
    fprintf(fid,'%s\t',title);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','Accuracy');
    fprintf(fid,'%g\t',acc);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','Precision');
    fprintf(fid,'%g\t',prec);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','Recall');
    fprintf(fid,'%g\t',rec);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','F1');
    fprintf(fid,'%g\t',f1);
    fprintf(fid,'\r\n'); % 换行

    fprintf(fid,'%s\t','AccStd');
    fprintf(fid,'%g\t',acc_std);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','PreStd');
    fprintf(fid,'%g\t',prec_std);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','RecStd');
    fprintf(fid,'%g\t',rec_std);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\t','F1Std');
    fprintf(fid,'%g\t',f1_std);
    fprintf(fid,'\r\n'); % 换行
    fprintf(fid,'%s\r\n',repmat('----', 1, 16)); % 分割线
    fclose(fid);
end