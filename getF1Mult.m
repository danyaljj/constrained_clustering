function [ f1arr ] = getF1Mult( X, Y )

k=length(unique(Y));

f1arr= []

%% cluster: dp-means: 
% lambda = 10; 
T = mean(X);
%[dist, ind] = sort( sqrt(sum((repmat(T,size(X,1),1)-X).^2,2)), 'descend' );
dists = []
for i=1:1:size(X,1)
    dists(end+1) = multDifference(T, X(i,:));
end
[dist, ind] = sort(dists, 'descend' );

lambda = dist(k);

[centroid, pointsInCluster, assignment, clusterSize]= dpmeans(X, lambda, 'Multinomial'); 

[p,r] = getbcubed(Y, assignment);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

end
