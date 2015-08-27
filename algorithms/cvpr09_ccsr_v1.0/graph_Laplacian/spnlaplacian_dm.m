function nL = spnlaplacian_dm(distMatrix,sigma,k)

% nL = nlaplacian_sdm(sdm,sigma)
% compute the normalized graph Laplacian from squared distance matrix
% distMatrix - squared distance matrix

nW = spnsimilarity_dm(distMatrix,sigma,k);

nL = speye(size(nW,1)) - nW;