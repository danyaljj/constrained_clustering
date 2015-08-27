function e = entroa(cmat)

% function e = entroa(cmat)
%
% Calculates the average entropy (e) of a clustering, given by the
% confusion matrix (cmat). Each column in the matrix corresponds
% a class and each row a cluster.

kc = size(cmat,1);
e0 = zeros(kc,1);
nc = sum(cmat,2);
I = find(nc>0);
for i = 1 : length(I)
    e0(I(i)) = entro(cmat(I(i),:));
end

e = mean(e0(I))/log2(size(cmat,2));

return;
