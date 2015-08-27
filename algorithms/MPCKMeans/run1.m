clear all; close all; clc;



% classifier = weka.clusterers.MPCKMeans;
% weka.clusterers.MPCKMeans.testCase()
% weka.clusterers.MPCKMeans.main('-D sdfsdf -C aeaqwe')
% rate = 0.01;
% [X,Y] = gaussians(200);
% size1 = size(X,1);
% order = randperm(size1);
% X = X(order,:);
% Y = Y(order,:);
% k = 2;    

% p = 1; 
% data = twospirals();
% X = data(:,1:2);
% Y = data(:,3);
% size1 = size(X,1);
% order = randperm(size1);
% X = X(order,:);
% Y = Y(order,:);
% k = 2;

k = 4;
%for k = 1:10
data = corners();
X = data(:,1:2);
Y = data(:,3);
size1 = size(X,1);
order = randperm(size1);
X = X(order,:);
Y = Y(order,:);

[assignment, centroids_vectors] = runMPCKMeans(X, Y, k, []); 


% for i=1:mpckmeans.getTotalTrainWithLabels.numInstances()
% %     i 
%     cluster_i = assignment(i);    
%     ins_i = totalTrainWithLabels.instance(i-1); % take care of java's 0-indexin pattern. 
%     class_i = ins_i.classValue();
%     for j=i+1:mpckmeans.getTotalTrainWithLabels.numInstances()
% %         j 
% %         tmp = mpckmeans.m_ClusterAssignments; 
%         cluster_j = assignment(j);
%         ins_j = totalTrainWithLabels.instance(j-1);  
%         class_j = ins_j.classValue();
%         if (cluster_i == cluster_j && class_i == class_j) || (cluster_i ~= cluster_j && class_i ~= class_j) 
%             nCorrect = nCorrect + 1;
%         end  
%     end 
% end 
% numInstances = mpckmeans.getTotalTrainWithLabels.numInstances();
% RandIndex = 100.0 * nCorrect/(numInstances*(numInstances-1)/2);
% disp(['Acc\t  '  num2str(RandIndex)]);
        

h = figure;
hold on;
% title(t)
colormap([1 0 .5;   % magenta
           0 0 .8;   % blue
           0 .6 0;   % dark green
           .3 1 0]); % bright green
dotsize = 12;
scatter(X(:,1), X(:,2), dotsize, assignment); axis equal;


