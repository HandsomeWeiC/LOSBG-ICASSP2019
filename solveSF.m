% min_{S>=0, S'*1=1, S*1=1, F'*F=I}  ||S - A||^2 + 2*lambda*trace(F'*Ln*F)
function [y,clusternum,S] = solveSF(A, c, islocal)
% Code for the following paper:
% Feiping Nie, Xiaoqian Wang, Cheng Deng, Heng Huang. 
% Learning A Structured Optimal Bipartite Graph for Co-Clustering.
% In the 31st Annual Conference on Neural Information Processing Systems (NIPS 2017).
%
% input:
%   A: n*d input bipartite graph
%   c: number of clusters
%
% output:
%   y1: cluster indicator of n samples
%   y2: cluster indicator of d features
%   S: learned bipartite graph with explicit cluster structure
%
% author: Feiping Nie

if nargin < 3
    islocal = 1;
end;

NITER = 30;
zr = 10e-5;
lambda = 1;

[n,m] = size(A);
d = sum(A,2); D = diag(d);
La = D-A;
[F,ev0,ev] = eig1(La,c,0);

F = F(:,1:c);

if sum(ev(1:c+1)) < zr
    error('The original graph has more than %d connected component', c);
end;

a(:,1) = ev;

idxa = cell(n,1);
for i=1:n
    if islocal == 1
        idxa0 = find(A(i,:)>0);
    else
        idxa0 = 1:m;
    end;
    idxa{i} = idxa0;
end;

idxam = cell(m,1);
for i=1:m
    if islocal == 1
        idxa0 = find(A(:,i)>0);
    else
        idxa0 = 1:n;
    end;
    idxam{i} = idxa0;
end;


for iter = 1:NITER
    dist = L2_distance_1(F',F');  % only local distances need to be computed. speed will be increased using C
    %S = sparse(n,m);
    %S = spalloc(n,m,10*5);
    S = zeros(n,n);
    for i=1:n
        idxa0 = idxa{i};
        ai = A(i,idxa0);
        di = dist(i,idxa0);
        ad = (ai-0.5*lambda*di); 
        S(i,idxa0) = EProjSimplex_new(ad);
    end;
    
    
    %Sm = sparse(m,n);
    Sm = zeros(m,n);
    for i=1:m
        idxa0 = idxam{i};
        ai = A(idxa0,i);
        di = dist(idxa0,i);
        ad = (ai-0.25*lambda*di); 
        Sm(i,idxa0) = EProjSimplex_new(ad);
    end;


    SS = (S+Sm')/2;
    d = sum(SS,2); D = diag(d);
    Ls = D-SS;
    [F,ev0,ev] = eig1(Ls,c,0);
    a(:,iter+1) = ev;
    F_old = F;
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.0000001
        lambda = 2*lambda;
    elseif fn2 < 0.0000001
        lambda = lambda/2;   F = F_old;
    else
        break;
    end;
end;

SS0=sparse(n+m,n+m); SS0(1:n,n+1:end)=SS; SS0(n+1:end,1:n)=SS';
[clusternum, y]=graphconncomp(SS0);
y1=y(1:n)';
y2=y(n+1:end)';
y = y2;
end



