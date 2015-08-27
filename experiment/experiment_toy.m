clear all; close all; clc;
pathAll(''); 

% This is the rate in which we sample side information; 
% more accurate, this is the size of the pairwise constraints 
% over the total number of constraints. 
rate = 0.01;

% The confidence on quality of the constraints. When p = 1, the constraints
% are high quality (no noise). As `p` gets closer to zero, it is more 
% likely that the constraints will be fliped (hence noisier). 
p = 1;

%% Mixture of Gaussians data 
[X,Y] = gaussians(200);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;    
runClustering(X,Y,k, 'Gaussian-Mixtures', rate, p);

%% two spirals data 
data = twospirals();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;
runClustering(X,Y,k, 'Two-Spirals', rate, p);

%% Cluster In Cluster dataset 
data = clusterincluster();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;
runClustering(X,Y,k, 'Cluster-In-Cluster', rate, p);

%% Corners dataset 
k = 4;
%for k = 1:10
data = corners();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
runClustering(X,Y,k, 'Corners', rate, p);

%% Half-kernels dataset 
data = halfkernel();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;
runClustering(X,Y,k, 'Half-Kernel', rate, p);

%% Full-moon dataset 
data = crescentfullmoon();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;
runClustering(X,Y,k, 'crescentfullmoon', rate, p);

%% Outlier 
k = 4;
data = outlier();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
runClustering(X,Y,k, 'outlier', rate, p);

