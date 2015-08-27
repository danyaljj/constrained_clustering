function [M C] = genPWC(labels,nM,nC)

% [M C] = genPWC(labels,nM,nC)
% generate pairwise constraints
% nM - number of must-link constraint to be generated for each cluster
% nC - number of cannot-link constraint to be generated for every two
% clusters

labels = labels(:);u = unique(labels);k = length(u);

M = []; C = [];

% must-link
for i = 1:k
    idx = find(labels==u(i));
    n = length(idx);
    M = [M;[idx(ceil(rand(nM,1)*n)),idx(ceil(rand(nM,1)*n))]];
end

% cannot-link
for i = 1:k
    for j = i+1:k
        idx1 = find(labels==u(i)); idx2 = find(labels==u(j));
        n1 = length(idx1); n2 = length(idx2);
        C =[C;[idx1(ceil(rand(nC,1)*n1)),idx2(ceil(rand(nC,1)*n2))]];
    end
end
