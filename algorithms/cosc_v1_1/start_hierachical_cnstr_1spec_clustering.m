function [best_cut, best_clusters, best_viols, best_clusters_intermediate] = ...
            start_hierachical_cnstr_1spec_clustering(W, deg, ML, CL, k, start_flags, nInitializationsForBinarySplit, nRunsOfRecursiveSplit, MAX_ITERS, verbosity)                 
% Performs Constrained 1-Spectral Clustering as described in the paper, 
%   Syama Sundar Rangapuram, Matthias Hein
%   Constrained 1-Spectral Clustering   
%   In International conference on Artificial Intelligence and Statistics (AISTATS), 2012. 
%   Available here: http://jmlr.csail.mit.edu/proceedings/papers/v22/sundar12/sundar12.pdf
%
% Usage:
%   [cut, clusters, viols] = cosc(W, deg, k, ML, CL, start_flags, nRuns, perturbation, MAX_ITERS, verbosity)
%
% Input
%   W:      n x n similarity matrix, where n is the number of data points.
%           It has to be a sparse symmetric matrix, with zeros on the diagonal.
%   deg:    n x 1 degree vector (can also be seen as vertex weights)
%           for Normalized cut it is sum(W,2)
%           for Ratio cut, it is ones(size(W,1),1)
%   ML:     m x 2 martix, specifying m must-link pairs
%   CL:     p x 2 matrix, specifying p cannot-link pairs
%   k:      number of clusters
% 
%   OPTIONAL arguments:
%   start_flags:    0, don't use any other special initializations 
%                   >=1, initialize the method with the second eigenvector of
%                      the standard Graph Laplacian (Default for k>2)
%                      We used this option in the paper for both k=2 and k>2
%                   >=2, initialize with the solution of the method: Fast Normalized
%                      Cuts with Linear Constraints (Default for k=2)
%                      Note that this option is valid only for the binary parititioning case, i.e., k=2)
%   nRuns:          how many total runs? We start the method from (nRuns - start_flags) random
%                   initializations.
%                   Recommendation: in the paper we used 10 (but more is better).
%   perturabation:  true/false. Perturbing the solution will lead to better
%                   results but at the expense of higher runtime. (Currently
%                   available only for binary  partitioning and by default true for k=2.)
%   MAX_ITERS:      Maximum number of FISTA iterations for solving the inner problem.
%                   Recommendation: 1000 - 10000, depending on the size of the problem.                 
%   verbosity:      Display Level
%                   0: Display nothing. (Default option)
%                   1: Track the progress by displaying the main steps.
%                   2: Display intermediate results (useful for preliminary sanity checks).
%                   3: Display all details. Use this option to check if the maximum number
%                      of FISTA iterations reached before the inner problem
%                      is solved optimally.
% 
% Output
%   cut:        cut value of the constrained cut problem
%   clusters:   Indicator vector giving the clustering of the data points 
%   viols:      number of violated constraints: zero for the binary
%               partitioning problem
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    nOuterRuns = nRunsOfRecursiveSplit;
    nInnerRuns = nInitializationsForBinarySplit;
    cuts = zeros(nOuterRuns,1);
    clusters = cell(nOuterRuns,1);
    viols = inf*ones(nOuterRuns,1);
    clusters_intermediate = cell(nOuterRuns,1);

    for i=1:nOuterRuns

        if verbosity >= 1
            fprintf('\n******************************************STARTING OUTER RUN %d******************************************\n', i);
        end
        
        [cuts(i), clusters{i}, viols(i), clusters_intermediate{i}] = ...
            hierarchical_cnstr_1spec_clustering(W,ML,CL,deg,k, start_flags, nInnerRuns-start_flags, MAX_ITERS, verbosity);

        if verbosity >= 1
            fprintf('The Multi-cut = %f\t\t Corresponding viols = %d\n', cuts(i), viols(i));
        end
%        [cut, clusters, viols, clusters_intermediate] = hierarchical_cnstr_1spec_clustering(W,ML,CL,deg,k, nInnerRuns, verbosity);     
    end

    [min_cut, min_cut_ix] = min(cuts);        
    best_cut = min_cut;
    best_clusters = clusters{min_cut_ix};
    best_viols = viols(min_cut_ix);
    best_clusters_intermediate = clusters_intermediate{min_cut_ix};
        
end