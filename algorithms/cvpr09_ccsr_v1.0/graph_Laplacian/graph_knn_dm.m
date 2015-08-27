function W = graph_knn_dm(dm,k)

% W = graph_knn_dm(dm,k)
% form symmetric knn graph
% dm - distance matrix
% k - number specifying the neighborhood 

Npts = size(dm,1);
W = spalloc(Npts, Npts, Npts * k);

for i = 1 : Npts
    [tmp,idx] = sort(dm(:,i));
    W(i,idx(2:k+1)) = 1;
end

W = W | W';