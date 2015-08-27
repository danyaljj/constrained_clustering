function [vmin2,cut2, clusters2, viols2] = ...
    solve_cnstr_functional_incremental(W,CL,vertex_weights,start_flags, nRuns, prev_clusters, perturbation, MAX_ITERS, verbosity)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    n = size(W,1);
    bs_mode = false;
    gammas = 0;
    if ~isempty(CL)
        Q = construct_cnstr_graph(CL, n);
    else
        Q = sparse(n,n);
    end

    starts = cell(nRuns+start_flags,1);
    vmins = cell(nRuns+start_flags, length(gammas));
    cuts=0; viols=0; objectives=0;

    if start_flags >=2 
        if isempty(CL)
            starts{2} = randn(n,1);
            start_flags = 1;
        else
            if verbosity >= 2
                fprintf('%2sInitialization: computing the solution of the method: Fast normalized cuts with linear constraints\n', ' ');
            end
            labels = spec_clustering_lcnstrs_merging(W, vertex_weights, CL, [], true, verbosity);
            %%% The above method also returns apart from the indicator vector, the "constrained eigenvector",
            %%% which equals to v2 (second eigenvector of the Laplacian) if there are no constraints. 
            starts{2} = labels;
        end
    end

    if start_flags >=1 
        if verbosity >= 2
            fprintf('%2sIntialization: computing the second eigenvector of the standard Laplacian\n',' ');
        end
        v2 = eig_std_Laplacian(W,true,vertex_weights);
        starts{1} = v2;
    end

    if nRuns > 0 && verbosity >= 2
        fprintf('%2sInitialization: choosing %d random points\n', ' ',nRuns);
    end

    for t=start_flags+1:nRuns+start_flags
        starts{t} = prev_clusters{t-start_flags};
    end

    j = 1;
    flag = true;
    fprintf('\n');
    if verbosity >= 1 && ~isempty(CL)
        fprintf('%2sIncreasing gamma until all the constraints are satisfied.\n', ' ');
    end
    
    while flag
        if verbosity >= 1
            fprintf('%2sgamma = %f\n', ' ', gammas(j));
        end
        for  t=1:nRuns+start_flags
            if verbosity >= 2
                fprintf('%4sInitialization: %d\n', ' ', t);
            end

            %cnstr_prob_start = tic;
            if j==1
                [vmins{t,j},fmins]=solve_cnstr_functional(W,Q,starts{t},vertex_weights,gammas(j),MAX_ITERS, verbosity);
            else
                [vmins{t,j},fmins]=solve_cnstr_functional(W,Q,vmins{t, j-1},vertex_weights,gammas(j),MAX_ITERS,verbosity);
            end
            %toc(cnstr_prob_start);

            clusters = opt_thresh_cnstr_functional(vmins{t,j}, W, vertex_weights', Q,gammas(j), 1);        
            all_clusters{t,j} = clusters;
            cuts(t,j) = bal_cut(W, vertex_weights, clusters);
    %         errs(t,j) = cluster_err(clusters, Y);
            
            viols(t,j) = 0;
            if ~isempty(CL)
                viols(t,j) = sum(clusters(CL(:,1)) == clusters(CL(:,2)));    
            end
            objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters);

            if verbosity >= 2
                fprintf('%6sObjective value (after optimal thresholding): %f\t', ' ', objectives(t,j));
                fprintf('\t Corresponding cut: %f\t Corresponding viols: %d\n', cuts(t,j), viols(t,j));
            end

            if perturbation
            % perturbation           
                start = vmins{t,j};
                pert_counter = 0;
                counter = 1;
                for p=0.001:0.01:0.10%p=0.1:0.1:0.3%p=0.001:0.01:0.1 %p=0.1:0.1:0.7
                    if verbosity >= 2
                        fprintf('%6sPerturbing the solution. Count: %d\n', ' ', counter);
                    end
                    counter = counter + 1;
                    %pstart = vmins{t,j} + p* range(vmins{t,j}) * randn(n,1);
                    pstart = start + p* range(start) * randn(n,1);
                    [vmin,fmins]=solve_cnstr_functional(W,Q,pstart,vertex_weights,gammas(j),MAX_ITERS,verbosity);
                    clusters = opt_thresh_cnstr_functional(vmin, W, vertex_weights', Q,gammas(j), 1);                                       

                    if fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters) < objectives(t,j)
                        
                        vmins{t,j} = vmin;
                        all_clusters{t,j} = clusters;
                        cuts(t,j) = bal_cut(W, vertex_weights, clusters);
        %                 errs(t,j) = cluster_err(clusters, Y);

                        viols(t,j) = 0;
                        if ~isempty(CL)
                            viols(t,j) = sum(clusters(CL(:,1)) == clusters(CL(:,2)));    
                        end
                        objectives(t,j) = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters);

                        if verbosity >= 3
                            fprintf('%6sPerturbation improves the result\n', ' ');
                        end
                        
                    else
                        pert_counter = pert_counter+1;
                    end
                    
                    if verbosity >= 2
                        vt = 0;
                        if ~isempty(CL)
                            vt = sum(clusters(CL(:,1)) == clusters(CL(:,2)));    
                        end
                        fprintf('%6sObjective value (after optimal thresholding): %f\t', ' ', fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gammas(j), clusters));
                        fprintf('\t Corresponding cut: %f\t Corresponding viols: %d\n', bal_cut(W, vertex_weights, clusters), vt);
                    end

                    if pert_counter >=5 
                        break;
                    end
                end
            end
        end

        [min_o, min_i] = min(objectives(:,j));
        if verbosity >= 1
            fprintf('\n%4sBest Result (in terms of the Objective value): \t Cut = %f\t', ' ', cuts(min_i,j));
            fprintf('\t Corresponding viols = %d\n\n', viols(min_i,j));
        end

        if viols(min_i,j) == 0
            vmin2 = vmins{min_i, j}; cut2 = cuts(min_i, j); clusters2 = all_clusters{min_i, j}; viols2 = viols(min_i,j);
            if length(gammas) == 1 || (gammas(end)-gammas(end-1))/gammas(end)<1e-2                
                break;
            else
                max_gamma = gammas(end);
                gammas(end) = 0.5*(gammas(end)+gammas(end-1));
                bs_mode = true;
                
            end
        else if j==length(gammas) % j is always same as length(gammas)
            if bs_mode
                if (gammas(end)-gammas(end-1))/gammas(end)<1e-2
                    break;
                else
                    gammas = [gammas 0.5*(gammas(end)+max_gamma)];
                end
            else
            %new_gammas = gammas(end)+1:scale:gammas(end)+10*scale;
            if gammas(end) > 0
                gammas = [gammas gammas(end)*2];
            else
                gammas = [gammas gammas(end)+0.2];
            end    
            end
            end
            j = j+1;

        end

        
    end

%      [min_o, min_i] = min(objectives(:,j));
%      vmin = vmins{min_i,j};
%      cut = cuts(min_i, j);
%      clusters = all_clusters{min_i, j};
%      viols = viols(min_i, j);
%    viols = 0;

end
