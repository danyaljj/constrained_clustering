function p = puritya(cmat)

% function p = puritya(cmat)
%
% Calculates the average purity (p) of a clustering, given by
% the confusion matrix (cmat). In the confusion matrix, we restrict
% that each column is a class and each row a cluster.

nc = sum(cmat,2);
mc = max(cmat,[],2);
I = find(nc>0);
mc(I) = mc(I) ./ nc(I);

p = mean(mc(I));

return;
