function val = averagekmin_dm(dm,k)

% val = averagekmin_dm(dm,k)
% compute the averaged distance of each point to its k-th nearest neigbor.
% dm - distance matrix of a data set

dm = sort(dm);

val = mean(dm(k+1,:));