function nW = spnsimilarity_dm(distMatrix,sigma,k)

% nW = spnsimilarity_dm(distMatrix,sigma,k)
% compute the sparse normalized similarity matrix from a distance matrix.
% X - data matrix whose rows corresond to data points
% sigma - scale factor 
% k - number of neighbors

W = spsimilarity_dm(distMatrix,sigma,k);

E = sum(W);
E = 1./sqrt(E);

nW = W.*(E'*E);