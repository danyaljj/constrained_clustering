clear all; close all; clc;

addpath('../distance')

%% generate sample continues data
size1 = 100; 
mean1 = [2,-1];
cov1 = [1,0.1; 0.1,1];
mean2 = [8,3];
cov2 = [1 .2; 0.2,1];
X = [mvnrnd(mean1, cov1, size1); mvnrnd(mean2, cov2, size1)];
Y = [ones(size1,1)  ; -1*ones(size1,1)]; 
order = randperm(2*size1); 
X = X(order,:); 
Y = Y(order,:); 
k = 2;

%% make E matrix
E = zeros(size(Y, 1), size(Y, 1)); 
Checked = zeros(size(Y, 1), size(Y, 1)); 
randSize = .01 * size(Y, 1) * size(Y, 1); 
iterAll = 1;
while(1)
    i1 = randi(size(Y, 1)); 
    i2 = randi(size(Y, 1)); 
    if i1 == i2
        continue;
    end
    
    if Checked(i1, i2) == 0   
        Checked(i1, i2) = 1;
        if Y(i1) == Y(i2)
            E(i1, i2) = 1; 
            E(i2, i1) = 1; 
        else 
            E(i1, i2) = -1; 
            E(i2, i1) = -1; 
        end
        iterAll = iterAll + 1;
    end 
    if( iterAll > randSize) 
        break;
    end
end

T = mean(X);
[dist, ind] = sort( sqrt(sum((repmat(T,size(X,1),1)-X).^2,2)), 'descend' );

lambda = dist(k);
[centroid, pointsInCluster, assignment_out, clustersSize, objs, pointsAll, centroindsAll] = constrained_v2_newInit(X, lambda, E, 2, .001, 'Gaussian');
%[centroid, pointsInCluster, assignment, clustersSize, objs, pointsAll, centroindsAll] = constrained_v2(X, lambda, E, 2, .001, 'Gaussian');

[centroid, pointsInCluster, assignment_out, clustersSize, objs, pointsAll, centroindsAll] = constrained_dpmeans_slow(X, lambda, E, 2, .001);



assignment_out

