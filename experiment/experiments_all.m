clear all; close all; clc;
pathAll(''); 
rate = 0.01;


for iii = 1:5
%for p = [0.8 0.9 0.95]
    p = 1; 
    [X,Y] = gaussians(200);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    k = 2;    
    experiment_1(X,Y,k, 'Gaussian-Mixtures', rate, p);

    p = 1; 
    data = twospirals();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    k = 2;
    experiment_1(X,Y,k, 'Two-Spirals', rate, p);
    
    data = clusterincluster();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    k = 2;
    experiment_1(X,Y,k, 'Cluster-In-Cluster', rate, p);
    
    k = 4;
    %for k = 1:10
    data = corners();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    %     k = 2;
    experiment_1(X,Y,k, 'Corners', rate, p);
    %end
    
    p = 1; 
    data = halfkernel();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    k = 2;
    experiment_1(X,Y,k, 'Half-Kernel', rate, p);
    
    data = crescentfullmoon();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    k = 2;
    experiment_1(X,Y,k, 'crescentfullmoon', rate, p);
    
    p=1; 
    k = 4;
    %for k = 1:10
    data = outlier();
    X = data(:,1:2);
    Y = data(:,3);
    size1 = size(X,1);
    order = randperm(size1);
    X = X(order,:);
    Y = Y(order,:);
    %k = ;
    experiment_1(X,Y,k, 'outlier', rate, p);
    %end
%end

end 
