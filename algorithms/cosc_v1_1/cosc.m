function [cut, clusters, viols, clusters_intermediate] = cosc(W, deg, k, ML, CL, start_flags, nInitializationsForBinarySplit, nRunsOfRecursiveSplit, perturbation, MAX_ITERS, verbosity)
% Performs Constrained 1-Spectral Clustering as described in the paper, 
%   Syama Sundar Rangapuram, Matthias Hein
%   Constrained 1-Spectral Clustering   
%   In International conference on Artificial Intelligence and Statistics (AISTATS), 2012. 
%   Available here: http://jmlr.csail.mit.edu/proceedings/papers/v22/sundar12/sundar12.pdf
%   Long version: http://jmlr.csail.mit.edu/proceedings/papers/v22/sundar12/sundar12Supple.pdf
%
% Usage:
%   [cut, clusters, viols, clusters_intermediate] = cosc(W, deg, k, ML, CL, start_flags, nInitializationsForBinarySplit, nRunsOfRecursiveSplit, perturbation, MAX_ITERS, verbosity)
%
% Input
%   W:      n x n similarity matrix, where n is the number of data points.
%           It has to be a sparse symmetric matrix, with zeros on the diagonal.
%   deg:    n x 1 degree vector (can also be seen as vertex weights)
%           for Normalized cut it is sum(W,2)
%           for Ratio cut, it is ones(size(W,1),1)
%   k:      number of clusters
%   ML:     m x 2 martix, specifying m must-link pairs
%   CL:     p x 2 matrix, specifying p cannot-link pairs
% 
%   OPTIONAL arguments:
%   start_flags:    0, don't use any other special initializations 
%                   >=1, initialize the method with the second eigenvector of
%                      the standard Graph Laplacian (Default for k>2)
%                      We used this option in the paper for both k=2 and k>2
%                   >=2, initialize with the solution of the method: Fast Normalized
%                      Cuts with Linear Constraints (Default for k=2)
%                      Note that this option is valid only for the binary parititioning case, i.e., k=2)
%   nInitializationsForBinarySplit: 
%                   how many total initializations for each binary split? We start the method from (nInitializationsForBinarySplit - start_flags) random
%                   initializations.
%                   Default: 10 for k=2, 2 for k>2.
%   nRunsOfRecursiveSplit: 
%                   Total number of re-runs of the complete recursive
%                   split. This is valid only for k>2.
%                   Default: 5 for k>2; for k=2, value passed in this argument is ignored.
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
%   clusters:   vector containing the clustering of the data points.
%   viols:      number of violated constraints: zero for the binary
%               partitioning problem
%   clusters_intermediate:
%               matrix of size n times k-1, containing in each column the intermediate clusters found in the recursive
%               split. For example, clusters_intermediate(:, end) is same
%               as clusters, clusters_intermedaite(:, 1) is the first
%               2-way split and clusters_intermediate(:, i) gives
%               clustering into i+1 components.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    assert(nargin >= 5, 'Input error. Usage: [cut, clusters, viols, clusters_intermediate] = cosc(W, deg, k, ML, CL)');
    assert(k<=size(W,1), 'Wrong usage. Number of clusters is larger than the size of the graph.');
    assert(isnumeric(W) && issparse(W),'Wrong usage. W should be sparse and numeric.');
    assert(sum(sum(W~=W'))==0,'Wrong usage. W should be symmetric.');
    assert(sum(diag(W))==0,'Wrong usage. Graph contains self loops. W has to have zero diagonal.');	
    assert(k>=2, 'Input error. The number of clusters should be at least 2');
    
    if exist('start_flags', 'var')
        assert(start_flags <= 2, 'Input error. start_flags must be smaller than or equal to 2');
        assert(start_flags >= 0, 'Input error. start_flags must be larger than or equal to 0');
        if exist('nInitializationsForBinarySplit', 'var')
            assert(nInitializationsForBinarySplit >= start_flags, 'Input error. nInitializationsForBinarySplit must be larger than or equal to start_flags');
            if k>2 && exist('nRunsOfRecursiveSplit', 'var')
                assert(nRunsOfRecursiveSplit >= 1, 'Input error. nRunsOfRecursiveSplit must be larger than or equal to 1 for k > 2');
            end
        end
    end    
    
    if ~exist('verbosity', 'var'), verbosity = 0; end
    if ~exist('MAX_ITERS', 'var'), if k == 2, MAX_ITERS = 10000; else MAX_ITERS = 1000; end, end
    if ~exist('perturbation', 'var'), if k == 2, perturbation = true; else perturbation = false; end, end
    if ~exist('nRunsOfRecursiveSplit', 'var'), nRunsOfRecursiveSplit = 5; end
    if ~exist('nInitializationsForBinarySplit', 'var'), if k == 2, nInitializationsForBinarySplit = 10; else nInitializationsForBinarySplit = 2; end, end
    if ~exist('start_flags', 'var'), if k == 2, start_flags = 2; else start_flags = 1; end, end        	
        
    if k==2
        display('Computing (constrained) bi-partitioning with the following options:');
        display(['start_flags=', num2str(start_flags), ' nInitializationsForBinarySplit=', num2str(nInitializationsForBinarySplit), ' nRunsOfRecursiveSplit=', num2str(1), ' perturbation=', num2str(perturbation), ' MAX_ITERS=', num2str(MAX_ITERS), ' verbosity=', num2str(verbosity)]);
        [cut, clusters, viols] = start_cnstr_1spec_clustering(W, deg, ML, CL, start_flags, nInitializationsForBinarySplit, perturbation, MAX_ITERS, verbosity);
        clusters_intermediate = clusters;
    else
        %[cut, clusters, viols, clusters_intermediate] = start_hierachical_cnstr_1spec_clustering(W, deg, ML, CL, k, start_flags, nRuns, verbosity);
        if start_flags > 1
            start_flags = 1;
            display('For k>2, start_flags should be either 0 or 1; resetting start_flags = 1');
        end
        display('Computing (constrained) multi-partitioning with the following options:');
        display(['start_flags=', num2str(start_flags), ' nInitializationsForBinarySplit=', num2str(nInitializationsForBinarySplit), ' nRunsOfRecursiveSplit=', num2str(nRunsOfRecursiveSplit), ' perturbation=', num2str(perturbation), ' MAX_ITERS=', num2str(MAX_ITERS), ' verbosity=', num2str(verbosity)]);
        [cut, clusters, viols, clusters_intermediate] = start_hierachical_cnstr_1spec_clustering(W, deg, ML, CL, k, start_flags, nInitializationsForBinarySplit, nRunsOfRecursiveSplit, MAX_ITERS, verbosity);
    end
   
    fprintf('\nFinal Result: Cut = %f\t\t Number of violated constraints = %d\n', bal_cut(W, deg, clusters), viols);
end