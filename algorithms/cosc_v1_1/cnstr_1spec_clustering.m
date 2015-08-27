function [cut, clusters, viols] = ...
    cnstr_1spec_clustering(W,ML,CL,vertex_weights,start_flags,nRuns,prev_clusters,perturbation,MAX_ITERS,verbosity)
% Performs Constrained 1-Spectral Clustering 
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    assert(isnumeric(W) && issparse(W),'Wrong usage. W should be sparse and numeric.');
    assert(sum(diag(W))==0,'Wrong usage. W should have zeros on the diagonal')
    
    cnstr1 = [];
    cnstr2 = [];
    if ~isempty(CL)
        cnstr1 = CL(:,1);
        cnstr2 = CL(:,2);
    end
    
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% MUST LINK CONSTRAINTS %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%
    % Handle must-link constraints via sparsification.

    if verbosity>=1 && ~isempty(ML), display('Processing Must links via sparsification'); end
    %display('Processing Must links via sparsification'); 
	[W, vertex_weights, map, prev_clusters] = process_mls( W, vertex_weights, ML, prev_clusters );
    cnstr1 = map(cnstr1);
    cnstr2 = map(cnstr2);             
    
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% CANNOT LINK CONSTRAINTS %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%    
    % Now the cannot link constraints
    % First derive must links from cannot links 
   
%     if verbosity>=2
%         display('Deriving Must links from cannot links'); 
%     end    
    
    [dML] = derive_mls_frm_cls( W, [cnstr1, cnstr2] );
    [W, vertex_weights, map_dml, prev_clusters] = merge( W, vertex_weights, dML, prev_clusters);
    W = triu(W);    % merging of edges might yield non-symmetric edge weights with a very small error (< 1e-15)
    W = W + W';
    cnstr1 = map_dml(cnstr1);
    cnstr2 = map_dml(cnstr2);     
    
    if size(W,1)==2     % Check if the constraints already specified the partition!
        clusters = [0; 1];
        cut = bal_cut( W, vertex_weights, clusters);

        % Recover the labels on the original graph via two maps.. first dml's and then ml's; 
        clusters = clusters(map_dml);
        clusters = clusters(map);
        
        viols = 0;
        
    else
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% FUNCTIONAL MINIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
        % Now the real work starts!                
        dCL = CL;
        % Remove redundant cannot link constraints!
        if ~isempty(cnstr1) > 0
           dCL = [cnstr1 cnstr2];
           dCL = sort(dCL, 2);
           dCL = unique(dCL, 'rows');
           if verbosity >=1, fprintf('Minimizing the Functional subject to Cannot link constraints using %d starting points\n', nRuns+start_flags); end
        else
            if verbosity >=1, fprintf('Minimizing the unconstrained Functional using %d starting points\n', nRuns+start_flags); end
        end
        
        [vmin, cut, clusters, viols] = ...
            solve_cnstr_functional_incremental(W,dCL,vertex_weights,start_flags, nRuns, prev_clusters, perturbation,MAX_ITERS, verbosity);        

        % Recover the labels on the original graph via two maps.. first dml's and then ml's; 
        clusters = clusters(map_dml);
        clusters = clusters(map);       

    end
    
end