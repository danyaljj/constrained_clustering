clear;
addpath(genpath('.'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load fisheriris
x=meas;
x = x(1:5, :); 
y=strcmp('setosa',species)*1 + strcmp('versicolor',species)*2 + strcmp('virginica',species)*3;
y = y(1:5, :); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=3; % number of cluster
option = struct('init',1,'alpha',1,'rho2',1000,'bal',0,'distance',1);
nbConst=20;
Xi=1;

% compute distances between objects
D=x*x';
N=diag(diag(D))*ones(size(D));
DistObj=sqrt(N+N'-2*D);

% Generate constraints
link=addConstraints(y,nbConst);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,BetP,J,ab]=CEVCLUS(DistObj,K,link,Xi);

%%%%%%%%%%%%%%%
rate = 0.1; 
p = 1; 
size1 = size(x,1);
order = randperm(size1);
k = 3;    
experiment_1(x,y,k, 'Gaussian-Mixtures', rate, p);


