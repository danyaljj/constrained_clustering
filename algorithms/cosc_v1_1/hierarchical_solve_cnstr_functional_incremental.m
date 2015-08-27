function [vmin,cut1, cut2, clusters, viol]= ...
    hierarchical_solve_cnstr_functional_incremental(W,CL,vertex_weights,start_flags, nRuns, ...
    min_ncl, max_ncl,MAX_RUNS,FISTA,MAX_ITERS,verbosity)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

n = size(W,1);
%gammas = 0:0.1:1;
if ~isempty(CL)
    Q = construct_cnstr_graph(CL, n);
else
    Q = sparse(n,n);
    % gammas = 0;
end

gammas = 0;
beta = 0.1;
start_j = 1;

starts = cell(nRuns+start_flags,1);
vmins = cell(nRuns+start_flags, length(gammas));
cuts=0; viols=0; objectives=0;

if start_flags >=2
    if verbosity >= 2
        fprintf('%4sInitialization: Computing the solution of the method: Fast normalized cuts with linear constraints\n', ' ');
    end
    % Note that we still call the version with merging to make it more
    % efficient. This is the only option even for multi-partitioning as
    % the derived must links are enforced anyway!
    labels = spec_clustering_lcnstrs_merging(W, vertex_weights, CL, [], true, verbosity);
    starts{2} = labels;
end

% if start_flags >=3
%     [clusters0,cuts0,cheegers0] = computeMultiPartitioning(W,true,2,true,10,1,1,false);
%     starts{3} = clusters0;
% end

if start_flags >=1
    EigTime = tic;
    if verbosity >= 2
        fprintf('%4sIntialization: Computing the second eigenvector of the standard Laplacian\n',' ');
    end
    %[v2, lambda2, one_C] = eig_std_Laplacian(W,true,vertex_weights);
    [v2, lambda2, one_C] = eigs_Laplacian( W, vertex_weights );
    EigTime = toc(EigTime);
    if isnan(sum(v2))
        EigTime = tic;
        if verbosity >= 2
            fprintf('%4sIntialization: Computing the second eigenvector of the standard Laplacian using matlab eigs (we converged to zero vector)\n',' ');
        end
        [v2, lambda2, one_C] = eig_std_Laplacian(W,true,vertex_weights);
        EigTime = toc(EigTime);
    end
    if verbosity >= 3
        display(['Time Taken for Eigenvector computation: ', num2str(EigTime)]);
        display(['Cut value of spectral clustering (without cannot-links): ', num2str(fctval_cnstr_one_spec_Q(W, vertex_weights, Q, 0, one_C, true))]);
    end
    starts{1} = v2;
end

%TODO: draw more random points choose best 5 among them!
if nRuns > 0 && verbosity >= 2
    fprintf('%4sInitialization: Choosing %d random points\n', ' ',nRuns);
end

for t=start_flags+1:nRuns+start_flags
    starts{t} = randn(n,1);
end


scale = 1.1;
j = 1;
min_flag = true;

if verbosity >= 1 && ~isempty(CL)
    fprintf('\n%4sIncreasing gamma until atmost %d constraints are violated or at most %d gamma runs.\n', ' ', min_ncl, MAX_RUNS);
end

while j<=MAX_RUNS
    if verbosity >= 1
        fprintf('%4sgamma = %f\n', ' ', gammas(j));
    end
    
    for  t=1:nRuns+start_flags
        if verbosity >= 2
            fprintf('%6sInitialization: %d\n', ' ', t);
        end
        
        InnerProblemTime = tic;
        if j==1
            [vmins{t,j},fmins]=hierarchical_solve_cnstr_functional(W,Q,starts{t},vertex_weights,gammas(j),FISTA,MAX_ITERS,verbosity);
        else
            [vmins{t,j},fmins]=hierarchical_solve_cnstr_functional(W,Q,vmins{t, j-1},vertex_weights,gammas(j),FISTA,MAX_ITERS,verbosity);
        end
        InnerProblemTime = toc(InnerProblemTime);
        if verbosity >= 3
            display(['Time Taken for one initialization: ', num2str(InnerProblemTime)]);
        end
        clusters = opt_thresh_cnstr_functional(vmins{t,j}, W, vertex_weights', Q,gammas(j), 1);
        all_clusters{t,j} = clusters;
        cuts(t,j) = bal_cut(W, vertex_weights, clusters);
        viols(t,j) = 0;
        if ~isempty(CL)
            viols(t,j) = sum(clusters(CL(:,1)) == clusters(CL(:,2)));
        end
        %objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), vmins{t,j});
        obj = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), vmins{t,j});
        objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters);
        %assert( abs(objectives(t,j)-obj) <1e-12 || objectives(t,j) <= obj);
        if verbosity >= 2
            fprintf('%6sObjective value (after optimal thresholding): %f\t', ' ', objectives(t,j));
            fprintf('\t Corresponding cut: %f\t Corresponding viols: %d\n', cuts(t,j), viols(t,j));
        end
        
        % perturbation
        %         for p=0.001:0.01:0.12%p=0.1:0.1:0.3%p=0.001:0.01:0.1 %p=0.1:0.1:0.7
        %             pstart = vmins{t,j} + p* range(vmins{t,j}) * randn(n,1);
        %             [vmin,fmins]=solve_cnstr_functional(W,Q,pstart,vertex_weights,gammas(j),true,false);
        %             clusters = opt_thresh_cnstr_functional(vmin, W, vertex_weights', Q,gammas(j), 1);
        %
        %             if fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters) < objectives(t,j)
        %
        %                 vmins{t,j} = vmin;
        %                 all_clusters{t,j} = clusters;
        %                 cuts(t,j) = bal_cut(W, vertex_weights, clusters);
        %                 viols(t,j) = 0;
        %                   if ~isempty(CL)
        %                     viols(t,j) = sum(clusters(CL(:,1)) == clusters(CL(:,2)));
        %                   end
        %                 %objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), vmins{t,j});
        %                 obj = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), vmins{t,j});
        %                 objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), vmins{t,j});
        %                 assert( abs(objectives(t,j)-obj) <1e-12 || objectives(t,j) <= obj);
        %             end
        %         end
    end
    
    % It is better to check cuts rather than viols to increase gamma value
    % as in our case the penalty term has additional balancing term...
    if j>1 && abs( min(cuts(:,j)) - min(cuts(:,j-1)) )  < 1e-8
        %gammas(j+1:end) = gammas(j+1:end) * scale;
        beta = beta*2;
    end
    
    [min_o, min_i] = min(objectives(:,j));
    if verbosity >= 1
        fprintf('\n%4sBest Result (in terms of the Objective value): \t Cut = %f\t', ' ', cuts(min_i,j));
        fprintf('\t Corresponding viols = %d\t target viols =  %d\n\n', viols(min_i,j), min_ncl);
        %         fprintf('\n current viols %d\t target viols %f and %f\n', viols(min_i,j),min_ncl,max_ncl);
        %         fprintf('\n current cut %f\n', cuts(min_i,j));
    end
    
    if min_flag && viols(min_i,j) <= min_ncl
        start_j = j;
        min_flag = false;
        %        break;
    end
    if viols(min_i,j) <= max_ncl
        break;
    else if j==length(gammas)
            %new_gammas = gammas(end)+1:scale:gammas(end)+10*scale;
            gammas = [gammas gammas(end)+beta];
        end
    end
    
    j = j+1;
    
end

% if j>MAX_RUNS % we dont satsify minimum required, so choose the unconstrained version
%     j = 1;
% end

cuts_cnstr = cuts(:, start_j:end);
min_cut = min(min(cuts_cnstr));
[ix, jx] = find(cuts_cnstr==min_cut);
min_i = ix(1); min_j = jx(1);
min_j = min_j + start_j-1;

% cuts_cnstr = cuts(:, start_j:end);
% [objs, ixs] = min(cuts_cnstr);
% [min_o, min_j] = min(objs);
% min_i = ixs(min_j);
% [min_o, min_i] = min(objectives(:,j));

vmin = vmins{min_i,min_j};
%fmins = objectives(min_i,j);
cut = cuts(min_i, min_j);
clusters = all_clusters{min_i, min_j};

cut1 = sum(sum(W(clusters==1, clusters~=1)))/sum(vertex_weights(clusters==1));
cut1 = full(cut1);
cut2 = sum(sum(W(clusters==1, clusters~=1)))/sum(vertex_weights(clusters~=1));
cut2 = full(cut2);

%assert( abs(cut - cut1 - cut2) < 1e-12);
viol = viols(min_i,min_j);
end
