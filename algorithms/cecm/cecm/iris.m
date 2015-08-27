clear;
addpath(genpath('.'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load fisheriris
x=meas;
y=strcmp('setosa',species)*1 + strcmp('versicolor',species)*2 + strcmp('virginica',species)*3;
n=length(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=3; % number of cluster
option = struct('init',1,'alpha',1,'rho2',1000,'bal',0,'distance',1);

nbConst=10;
noise=0;
matConst=eye(n);
matConst=addNewConstraints(x,y,matConst,nbConst,noise,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,g,BetP,J]=CECM(x,K,matConst,option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 1; 
% [X,Y] = gaussians(200);
[X,Y] = twospirals();
% [X,Y] = clusterincluster();
    
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);
k = 2;
n=length(Y);
matConst=eye(n);
matConst=addNewConstraints(X,Y,matConst,nbConst,noise,0);
[m,g,BetP,J]=CECM(X,k,matConst,option);
[C,Idx] = max(BetP,[],2); 

scatter(X(:,1), X(:,2), 12, Idx); axis equal;
    