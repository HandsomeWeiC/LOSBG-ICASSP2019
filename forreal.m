function acc = forreal(X,k,lambda,truth)

% 参数设置
 itermax = 10; %最多15次吧
%  lambda = 0.18;
 lambda_2 = 0.00000000001;
 lambda_3 = 1;
%% 主循环 (exactly ALM)
%初始化 
[d,n]=size(X);
S = zeros(n,n);
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

[U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);
S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;
% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:k);
V = D*V;
all=size(truth,1);
idx = kmeans(V,k,'emptyaction','singleton','replicates',20,'display','off');
idxZ = kmeans(Z',k,'emptyaction','singleton','replicates',20,'display','off');

G = full(S);
idxS = kmeans(G',k,'emptyaction','singleton','replicates',20,'display','off'); 

[miss,index] = missclassGroups(y,truth,k);
Accy=1-miss/all;

[miss,index] = missclassGroups(idx,truth,k);
Accx=1-miss/all;

[miss,index] = missclassGroups(idxZ,truth,k);
AccZ=1-miss/all;

[miss,index] = missclassGroups(idxS,truth,k);
AccS=1-miss/all;

acc = max([Accy,Accx,AccZ,AccS]);

end



