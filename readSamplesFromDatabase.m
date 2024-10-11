function [sample,sample1,realclass] = readSamplesFromDatabase(databaseType,readNames,noise_rate,noise_type,imageRow,imageColumn)
close all;
%该函数用来从数据库中读取样本，读取的样本图片都为32 x 32维。

%输入参数:databaseType 表示要读取的数据库类型。
% readNames是要读取的每类样本的名字的矩阵，如该矩阵为 [7 10]，则表示要读取每一类的第七张和第十张图片。

%返回为样本矩阵。该矩阵维度为 总的读取样本数 x 每个样本维度(32 * 32)
% folder_now = 'E:\code for PCA';  addpath([folder_now, '\funs']);
if databaseType == 1 %3 5 9 11  k=3
	myDataPath = ('E:\code for PCA\yale\s');
    typeName = '.bmp';
	classCounts = 15;
elseif databaseType == 2
	typeName = '.jpg';
	classCounts = 100;
    dataPath = ['E:\code for PCA\AR\p'];
elseif databaseType == 3 
    typeName = '.pgm';
    classCounts = 38;
    dataPath = ['E:\code for PCA\CroppedYale\yaleB'];
elseif databaseType == 4   %k=4
    typeName = '.pgm';
    classCounts = 15;
    dataPath = ['E:\code for PCA\Umist\'];
elseif databaseType == 5   %20
    typeName = '.tiff';
    classCounts = 10;
    myDataPath = ['E:\code for PCA\Jaffe\s'];
elseif databaseType == 6
    typeName = '.jpg';
    classCounts = 20;
    dataPath = ['E:\code for PCA\Coil20\s'];
elseif databaseType == 7
    typeName = '.jpg';
    classCounts = 10;
    dataPath = ['E:\code for PCA\CBCL\s'];
elseif databaseType == 8
    typeName = '.jpg';
    classCounts = 50;
    myDataPath = ['E:\code for PCA\FEI\'];
elseif databaseType == 9
	typeName = '.jpg';
	classCounts = 100;
    dataPath = ['E:\code for PCA\AR_Data\AR_Data\AR\s'];
elseif databaseType == 10
	typeName = '.bmp';
	classCounts = 40;
    myDataPath = ['E:\code for PCA\ORL\s'];
elseif databaseType == 11
	typeName = '.jpg';
	classCounts = 68;
    dataPath = ['E:\code for PCA\CMU PIE pose09\s'];
    elseif databaseType == 12
	typeName = '.tif';
	classCounts = 200;
    dataPath = ['E:\code for PCA\FERET\FERET'];
    elseif databaseType == 13
	typeName = '.jpg';
	classCounts = 15;
    dataPath = ['E:\code for PCA\15_s\s'];
    elseif databaseType == 14
	typeName = '.jpg';
	classCounts = 1;
    dataPath = ['E:\code for PCA\indian pines\'];
    elseif databaseType == 15
	typeName = '.jpg';
	classCounts = 50;
    myDataPath = ['E:\code for PCA\GTdb_crop\cropped_faces\s'];
   elseif databaseType == 16
    typeName = '.jpg';
	classCounts = 20;
    dataPath = ['E:\code for PCA\pubfig83\s'];
    elseif databaseType == 17
    typeName = '.png';
	classCounts = 100;
    dataPath = ['E:\code for PCA\coil100\s'];
    elseif databaseType == 18
    typeName = '.jpg';
	classCounts = 20;
    dataPath = ['E:\code for PCA\101_OC\s'];
end



sample = zeros(classCounts * size(readNames,2), imageRow * imageColumn);
realclass = zeros(classCounts * size(readNames,2),1);
d = floor((classCounts * size(readNames,2))*noise_rate);%噪声比例
noise_index = randperm((classCounts * size(readNames,2)),d);

for i = 1:classCounts
    count = 0;
    
    if databaseType == 2 || databaseType == 3 || databaseType == 4||...
            databaseType == 9|| databaseType == 6||....
            databaseType == 7||databaseType == 11 ||databaseType == 16||databaseType == 17....
        ||databaseType == 18
        myDataPath = strcat(dataPath,num2str(i),'\');
    elseif databaseType == 12
          myDataPath = strcat(dataPath,'-',num2str(i),'\'); 
       elseif databaseType == 13
          myDataPath = strcat(dataPath,num2str(i),'\');
    end
    
    index = 1;
    
%       for ii=1:7
    for readNameNumber = readNames
        
        %path为读取的图片路径
        if databaseType == 1
            path = strcat(myDataPath,num2str(i),'_',num2str(readNameNumber),typeName);
        elseif databaseType == 2
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 3
            path = strcat(myDataPath,'(',num2str(readNameNumber),')',typeName);
        elseif databaseType == 4
            path = strcat(myDataPath,'(',num2str(readNameNumber),')',typeName);
        elseif databaseType == 5
            path = strcat(myDataPath,num2str(i),'_ (',num2str(readNameNumber),')',typeName);
        elseif databaseType == 6
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 7
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 8
             if readNameNumber < 10
                readNameNumberString = strcat('0',num2str(readNameNumber));
            else
                readNameNumberString = num2str(readNameNumber);
            end
            path = strcat(myDataPath,num2str(i),'-',num2str(readNameNumberString),typeName);
           
        elseif databaseType == 9
             path = strcat(myDataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 10
            path = strcat(myDataPath,num2str(i),'_' ,num2str(readNameNumber),typeName);
        elseif databaseType == 11
             path = strcat(myDataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 12
                readNameNumberString = strcat('0',num2str(readNameNumber));
                path = strcat(myDataPath,num2str(readNameNumberString),typeName);
        elseif databaseType == 13
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
             elseif databaseType == 14
            path = strcat(dataPath,num2str(readNameNumber),typeName);
        elseif databaseType == 15
            if readNameNumber < 10
                readNameNumberString = strcat('0',num2str(readNameNumber));
            else
                readNameNumberString = num2str(readNameNumber);
            end
            if i < 10
                iiii = strcat('0',num2str(i));
            else
                iiii = num2str(i);
            end
            path = strcat(myDataPath,num2str(iiii),'_',num2str(readNameNumberString),typeName);
            elseif databaseType == 16
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
             elseif databaseType == 17
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
            elseif databaseType == 18
            path = strcat(myDataPath,num2str(readNameNumber),typeName);
        end
        
         readData = imread(path);
%          readData= histeq(readData,255);
%           m = size(readData,3);
%           if m >1
%          readData = rgb2gray(readData);
%           end
%         imshow(readData);
%       imshow(path);
        readData = imresize(readData,[imageRow,imageColumn]);     
        readData1 = imresize(readData,[imageRow,imageColumn]);
% % %         max(max(readData))
% % %         min(min(readData)) 
       %subplot(222);  
        %imshow(readData1);
       %title('加入高斯噪声后的图像'); 
% % %         subplot(223);  
% % %         imshow(readData2);
% % %         title('加入椒盐噪声后的图像'); 
% % %         subplot(224);  
% % %         imshow(readData3);
% % %         title('加入乘性噪声后的图像'); 
        %%
%           subplot(212); 

        imageData1 = readData1(1:imageRow * imageColumn);
        line = (i-1) * size(readNames,2) + index;
        sample1(line,:) = imageData1; 

        h = floor(d*0.5);
        one_noise = noise_index(:,1:h);
        two_noise = noise_index(:,(h+1):d);
        if find(line==one_noise)
         [readData] = addnoise(readData,noise_type);
         imageData = readData(1:imageRow * imageColumn);
        elseif find(line==two_noise)
            [readData] = addnoise(readData,noise_type);
         imageData = readData(1:imageRow * imageColumn);
        else
         imageData = readData(1:imageRow * imageColumn);
        end
        imshow(readData);
        line = (i-1) * size(readNames,2) + index;
        sample(line,:) = imageData; 
        realclass(line,:) = i;
        index = index + 1;
        count = count + 1;
    end
end
end