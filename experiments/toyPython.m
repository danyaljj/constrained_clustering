clear all; close all; clc;

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

%% 
system('python var_dpmm.py')

