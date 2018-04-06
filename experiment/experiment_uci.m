% this script runs clustering on UCI datasets
clear all; close all; clc;
pathAll('');
% This is the rate in which we sample side information;
% more accurate, this is the size of the pairwise constraints
% over the total number of constraints.
rate = 0.03;

% The confidence on quality of the constraints. When p = 1, the constraints
% are high quality (no noise). As `p` gets closer to zero, it is more
% likely that the constraints will be fliped (hence noisier).
p = 1;


% rand of values for k 
k_variance=1

% number of trials
iii=3

%% iris 
[X,Y] = readUCIData('iris');
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 3;
experimentName = ['Variable_k_iris_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance) ];
runClustering(X,Y,k, experimentName, rate, p, false);

%% wine
[X,Y] = readUCIData('wine');
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 3;
experimentName = ['Variable_k_wine_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)];
runClustering(X,Y,k, experimentName, rate, p, false);

%% ecoli
[X,Y] = readUCIData('ecoli');
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 8;
experimentName = ['Variable_k_eclli_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)];
runClustering(X,Y,k, experimentName, rate, p, false);

%% glass
[X,Y] = readUCIData('glass');
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 7;
experimentName = ['Variable_k_glass_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)];
runClustering(X,Y,k, experimentName, rate, p, false);

%% balance
[X,Y] = readUCIData('balance');
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 3;
experimentName = ['Variable_k_balance_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)];
runClustering(X,Y,k, experimentName, rate, p, false);
