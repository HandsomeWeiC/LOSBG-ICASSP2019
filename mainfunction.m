% the optimal problem
% min_{Z,E,S,F} ||Z||_*+ \lambda_1||E||_{2,1}+\lambda_2||S-Z||_F^2...
% + \lambda_3 tr(F'L_SF)
% s.t. X=XZ+E, S>=0, S'1=1, F'F=I, F \in R^{n*k}.
% X: the data matrix;  Z: the coefficient matrix; E: the error matrix
% S; the approximat matrix to Z; L_S: the Laplacian matrix of S
tic
clear
clc
addpath('data')
%% 数据输入
X = cell2mat(struct2cell(load('X_highdimen0.mat'))); % data matrix
k = 5; % the cluster number
truth = cell2mat(struct2cell(load('truth200.mat')));
% 参数设置
 itermax = 15; %max
 lambda = 0.18;
 lambda_2 = 0.000001;
 lambda_3 = 1;
%% 主循环 (exactly ALM)
%初始化 
[d,n]=size(X);
S = zeros(n,n);
% mian function: LOSBG
 for iter=1:itermax
    % update Z,E
    [Z,E] = solve_lrr(X,S,lambda,lambda_2);

    % update S,F
%     [y,clusternum,S] = solveSF(Z,k);
    [y,clusternum,S] = coclustering_bipartite_fast(Z,k,lambda_3);
    if clusternum == k
        break;
    end

 end
 
%% 数据输出
[miss,index] = missclassGroups(y,truth,k);
all=size(truth,1);
Accy=1-miss/all
G = full(S);
% imagesc(G);
clean(G,1);

% imagesc(Z)
% clean(Z,2);

if clusternum ~= k
    sprintf('Can not find the correct cluster number: %d', k);
    idxZ = kmeans(Z',k,'emptyaction','singleton','replicates',20,'display','off');
    idxS = kmeans(G',k,'emptyaction','singleton','replicates',20,'display','off'); 
    [miss,index] = missclassGroups(idxZ,truth,k);
    AccZ=1-miss/all
    [miss,index] = missclassGroups(idxS,truth,k);
    AccS=1-miss/all
end;
toc


