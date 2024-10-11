
% EigenGap = load('coil20_convergence.mat');

%{
EigenGap= jjj;
plot(EigenGap,'-o','color','[0 0.7 0]','LineWidth',2,'MarkerEdgeColor','[0 0.7 0]','MarkerFaceColor','[0 0.7 0]');
xlabel('Iterations'),ylabel('EigenGap(c)')
axis([0 30 0 max(EigenGap)+0.1*max(EigenGap)])
%}
clear;clc;close all;

load Cars.mat

num = size(X, 1);
c = length(unique(Y));
feature_num = c-1;
k1 = floor(0.8*num);
beta = 1;
k = 5;
[W,obj,S,p, fn]=FS_AN(X',c,feature_num,k1,beta,k);
%obj = obj(2:30);
plotObj = zeros(50, 1);
for i = 1:50
    if obj(i) ~= 0
        plotObj(i) = obj(i);
    else
        plotObj(i) = plotObj(i-1);
    end
end
%%
plot(plotObj,'-o','color','r','LineWidth',2,'MarkerEdgeColor','r');
xlabel('Iterations', 'Fontsize', 24)
ylabel('The Objective Function', 'Fontsize', 24)
ax = gca;
ax.FontSize = 16;