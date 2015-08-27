function W = spsimilarity_dm(distMatrix,sigma,k)

% W = similarity_dm(dm,sigma)
% compute the sparse similarity matrix from a distance matrix using Gaussian kernel 
% w_ij = exp{-d^2(x_i,x_j)/(2*sigma^2)}
% distMatrix - a distance matrix
% sigma - scale factor
% k - number of neighbors

W = exp(-(distMatrix.^2)/(2*sigma^2));
A = graph_knn_dm(distMatrix,k); % adjacency matrix
W = A .* W;

Npts = size(distMatrix,1);
W(1:Npts+1:end) = 0;


