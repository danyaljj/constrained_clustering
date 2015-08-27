load two_moons.mat % this has the weight matrix W, must-links ML, cannot-links CL and ground truth Y (used only for computing the clustering error)

vertex_weights = sum(W,2); % this choice corresponds to normalized cut
%vertex_weights = ones(size(W,1),1); % this choice corresponds to ratio cut
k = length(unique(Y)); % no. of clusters


%------------------------------------------- no. of clusters = 2 ---------------------------------------------------%
% This is the default call. It starts the method from 10 different
% initializations and takes the best (according to the cut value) among them.
[cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL);  

% This is a fast version (with some compromise on quality). It starts the method from only 2 special
% initializations.
% [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000);  

% Use this to track progress by displaying the intermediate results
% [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000, 2);  

% To get the detailed help, type at the matlab command prompt: help cosc.m

%--------------------------------------------------------------------------------------------------------------------%


%------------------------------------------- no. of clusters > 2 ---------------------------------------------------%
% k = 4; 
% This is the default call. It takes the best result (according to the
% multi-cut) out of 5 complete recursive splits. In each step, it uses two initializations for
% computing the 2-way split.
% [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL);  

% This is a fast version (with some compromise on quality). It computes only one complete recursive split and in each step only 1 special
% initialization is used for the 2-way split.
% [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 1, 1, 1, false, 1000);  

% Use this to obtain the intermediate clusterings of the recursive split in the last output argument.
% For example, clusters_intermediate(:, i) provides the clustering result for i+1 clusters.
% [cut, clusters, viols, clusters_intermediate] = cosc(W, vertex_weights, k, ML, CL, 1, 1, 1, false, 1000);  

% Use this to track progress by displaying the intermediate results
% [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 1, 1, 1, false, 1000, 2);  

% To get the detailed help, type at the matlab command prompt: help cosc.m

%--------------------------------------------------------------------------------------------------------------------%
