function [cut, clusters, viols, clusters_intermediate_orig] = ...
    hierarchical_cnstr_1spec_clustering(W,ML,CL,vertex_weights,k, start_flags, nRuns, MAX_ITERS,verbosity)
% Performs Hierarchical Constrained 1-Spectral Clustering yielding
% multi-partitioning of the data.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FISTA = true;
MAX_RUNS = 5;
%start_flags = 0;
totCnstrs = size(ML,1) + size(CL,1);
n = size(W,1);
clusters_intermediate_orig=zeros(n,k-1);

cnstr1 = [];
cnstr2 = [];
if ~isempty(CL)
    cnstr1 = CL(:,1);
    cnstr2 = CL(:,2);
end

%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% MUST LINK CONSTRAINTS %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%
% Handle the must-link constraints via sparsification.
%ml_start = tic;
prev_clusters = cell(1,1);
prev_clusters{1,1} = rand(n,1);
%map = 1:size(W,1);
if verbosity >= 1 && ~isempty(ML), display('Processing Must links via sparsification'); end
MLTime = tic;
[W, vertex_weights, map, prev_clusters] = process_mls( W, vertex_weights, ML, prev_clusters );
MLTime = toc(MLTime);
if verbosity >= 1
    display(['Time Taken for Must links processing: ', num2str(MLTime)]);
end
cnstr1 = map(cnstr1);
cnstr2 = map(cnstr2);
W = triu(W);
W = W + W';
%toc(ml_start);

%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% CANNOT LINK CONSTRAINTS %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%
% Now the cannot link constraints
if size(W,1)==k     % Check if the constraints already specified the partition!
    clusters = 1:k;
    cut = bal_cut( W, vertex_weights, clusters);
    clusters = clusters(map);
    clusters_intermediate_orig = clusters;
    
    nVoilated = 0;
    if ~isempty(CL)
        nVoilated = sum( clusters(CL(:,1)) == clusters(CL(:,2)) );    % get the number violated on the original graph.
    end
    
    viols = nVoilated;
else
    %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% FUNCTIONAL MINIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
    % Now the real work starts!
    n = size(W,1);
    % Remove redundant cannot link constraints!
    if ~isempty(cnstr1)
        CL = [cnstr1 cnstr2];
        CL = sort(CL, 2);
        CL = unique(CL, 'rows');
        
        % Little sanity check to see we did not remove non-redundant
        % cannot links!
        %q1 = construct_cnstr_graph([cnstr1 cnstr2], n);
        %q2 = construct_cnstr_graph(CL, n);
        %assert( sum(sum(q1-q2)) == 0 );
    end
    
    %cnstr_prob_start = tic;
    knew = k;
    %min_ncl = size(CL,1)/3;
    min_ncl = 0.5*(knew - 2)*size(CL,1)/(knew-1);
    min_ncl = floor(min_ncl);
    min_ncl = (knew - 2)*size(CL,1)/knew;
    min_ncl = floor(min_ncl);
    k1 = floor(knew/2);
    k2 = knew-k1;
    max_ncl = size(CL,1)* (k1*(k1-1)+k2*(k2-1))/(knew *(knew-1));
    max_ncl = floor(max_ncl);
    max_ncl = min_ncl;
    if verbosity >=1
        if ~isempty(CL), fprintf('%2sMinimizing the Functional subject to Cannot link constraints using %d starting points\n', ' ', nRuns+start_flags);
        else fprintf('%2sMinimizing the unconstrained Functional using %d starting points\n', ' ', nRuns+start_flags);
        end
    end
    %min_ncl = size(CL,1);
    %[vmin,cut, clusters, vmins, cuts_all, objectives, viols]= ...
    %[vmin, cut, clusters, viols] = hierarchical_solve_cnstr_functional_incremental(W,CL,vertex_weights,start_flags, nRuns, min_ncl, FISTA,verbosity>=3);
    [vmin, cut1, cut2, clusters, viols] = hierarchical_solve_cnstr_functional_incremental(W,CL,vertex_weights,start_flags, nRuns, ...
        min_ncl, max_ncl,MAX_RUNS, FISTA,MAX_ITERS,verbosity);
    %toc(cnstr_prob_start);
    
    cut = cut1+cut2;
    mcut = bal_cut( W, vertex_weights, clusters);
    if cut ~= inf
        %assert( abs(mcut - cut) < 1e-12);
    end
    
    nsatis_tot = 0;
    if ~isempty(CL)
        nsatis_tot = sum(clusters(CL(:,1)) ~= clusters(CL(:,2)));
    end
    
    ncl = size(CL,1);
    
    %         if verbosity >= 3
    %             fprintf('Best result:\n');
    %             fprintf('Normalized Cut: %.8g \n\n',cut);
    %         end
    
    deg = sum(W,2);% do we need it???
    n = size(W,1);
    Q = construct_cnstr_graph(CL, n);
    
    clusters_intermediate=zeros(n,k-1);
    cuts_intermediate=zeros(1,k-1);
    cutParts=zeros(1,k);
    
    clusters=clusters+1;
    clusters_intermediate(:,1)=clusters;
    
    cuts_intermediate(:,1)=cut;
    
    cutPart1 = sum(sum(W(clusters==1, clusters~=1)))/sum(vertex_weights(clusters==1));
    cutPart2 = sum(sum(W(clusters==1, clusters~=1)))/sum(vertex_weights(clusters~=1));
    
    %assert( abs(cut - cutPart1 - cutPart2) < 1e-12 );
    cutParts(1)=cutPart1;
    cutParts(2)=cutPart2;
    
    sub_cuts=zeros(k,2);
    sub_clusters=cell(1,k);
    sub_nsatis=zeros(k,1);
    sub_nCL = zeros(k,1);
    
    if(verbosity>=1)
        fprintf('%2sFinished Clustering into 2 parts.  Sizes: %d, %d\n', ' ', sum(clusters==1), sum(clusters==2));
        fprintf('\n%2sStarting hierarchical parititioning...\n', ' ');
    end
    
    connectedComponentCounter = 0;
    result_conn_comp_vmin = cell(0);
    result_ix_conn_comp = sparse(size(W,1), size(W,1));
    % search for index_m and replace by ix_m
    % replace connectedcomponents function by our version
    % in each if we are updating cuts after dividing m and clusters
    % we need cutpart1 and cutpart2 from the first split!
    % when do we get disconnected components and when does deg(ix(m1))
    % is zero!
    for l=3:k
        best_cut = inf;
        best_satis = 0;
        current_best_nsatis = 0;
        current_best_m =1;
        %             if l==k-1 && verbosity >=3
        %                 disp('stop');
        %             end
        % in each step consider each of the current l-1 clusters and
        % divide the one the yields good cut value
        for m=1:l-1
            
            ix_m=find(clusters==m);
            
            % if we have already solved this subproblem
            if (~isempty(sub_clusters{m}))
                clusters_m = sub_clusters{m};
                cut_m1 = sub_cuts(m,1);
                cut_m2 = sub_cuts(m,2);
                nsatis = sub_nsatis(m);
                nCL_m = sub_nCL(m);
                if verbosity >= 2
                    display(['Subgraph size: ', num2str(length(clusters_m))]);
                end
                
                % if the current cluster has size 1 it cannot be further divided
            elseif(length(ix_m)==1)
                clusters_m=[];
                cut_m1 = inf;
                cut_m2 = inf;
                sub_clusters{m}=clusters_m;
                sub_cuts(m,1) = cut_m1;
                sub_cuts(m,2) = cut_m2;
                nsatis = 0;
                
                % if the current cluster has just two vertices!
            elseif(length(ix_m)==2)
                clusters_m=[0;1];
                sub_clusters{m}=clusters_m;
                
                cut_m1 = deg(ix_m(1))-W(ix_m(1),ix_m(1));
                cut_m2 = deg(ix_m(2))-W(ix_m(2),ix_m(2));
                
                if cut_m1 > 0
                    cut_m1 = cut_m1/vertex_weights(ix_m(1));
                end
                
                if cut_m2>0
                    cut_m2=cut_m2/vertex_weights(ix_m(2));
                end
                
                %cut_m = cut_m1+cut_m2;
                sub_cuts(m,1) = cut_m1;
                sub_cuts(m,2) = cut_m2;
                sub_nsatis(m) = 0;
                
                nsatis = 0;
                if ~isempty( find( (CL(:,1)==ix_m(1) & CL(:,2)==ix_m(2)) | (CL(:,1)==ix_m(2) & CL(:,2)==ix_m(1)), 1 ) )
                    nsatis = 1;
                end
                                
                % otherwise we have to compute the partition
            elseif length(ix_m)>2
                if(verbosity >=1 ), fprintf('\n%2sComputing partitioning of subgraph %d.\n',' ' ,m);end
                
                % extract subgraph and its connected components
                Wm=W(ix_m,ix_m);
                vertex_weights_m = vertex_weights(ix_m);
                size_m=size(Wm,1);
%                 if (size_m == 0)
%                     display('bad component size');
%                 end
                [comp,connected,sizes]=connectedComponents(Wm);
                if verbosity >= 2
                    display(['Number of connected components of this partition: ', num2str(length(sizes))]);
                end
%                 if size(sizes,1) == 1
%                     display(['sizes: ', num2str(sizes)]);
%                 else
%                     display(['sizes: ', num2str(sizes')]);
%                 end
                Q_m = Q(ix_m, ix_m);
                [cix, cjx] = find(triu(Q_m));
                CL_m = [cix cjx];
                
                if(verbosity>=3 && ~connected), fprintf('%2sSubgraph is not connected.\n',' '); end
                cut_m1 = inf;
                cut_m2 = inf;
                
                % if subgraph is not connected
                if (~connected)
                    
                    %check partition which has connected component as one
                    %cluster and rest as second cluster
                    try
                        for m1=1:length(sizes)
                            %if(true)
                            if(cut_m1+cut_m2 > 0)
                                if(verbosity>=2) fprintf('...Checking partition found by isolating connected component %d of %d.\n',m1,length(sizes)); end
                                
                                clusters_m_temp = double(comp==m1); % clustering of the entire subgraph that has m connected components: one on the m1 component, zero everwhere else on the subgraph
                                
                                clusters_m2 = zeros(size(clusters,1),1); % to compute the cut of the m1 component to the rest of the graph, we form the clustering of the whole graph with 1 on the m1 component, 0 everywhere else on the graph
                                clusters_m2(ix_m )= clusters_m_temp;
                                %cut_m2_temp = bal_cut(W, vertex_weights, clusters_m2); %computeCutValue(cluster_m2,W,normalized,deg2);
                                cut_m2_temp = sum(sum(W(clusters_m2==0,clusters_m2==1)))/sum(vertex_weights(clusters_m2==1));
                                
                                clusters_m1=zeros(size(clusters,1),1); % Now we compute the cut of the remaining components (all of them as one cluster) of the subgraph: note that the split is m1 vs rest of the components
                                clusters_m1(ix_m)=(clusters_m_temp==0);
                                %cut_m1_temp = bal_cut(W, vertex_weights, clusters_m1); %computeCutValue(cluster_m1,W,normalized,deg2);
                                cut_m1_temp = sum(sum(W(clusters_m1==0,clusters_m1==1)))/sum(vertex_weights(clusters_m1==1));
                                
                                if verbosity >= 2
                                    display(['cut part1: ', num2str(cut_m1_temp), ' cut part2: ', num2str(cut_m2_temp)]);
                                end
                                cut_temp = sum(cutParts)-cutParts(m)+cut_m1_temp+cut_m2_temp;
                                % Display current objective
                                if(verbosity>=2)
                                    displayCurrentObjective(cut_temp,-1,true);
                                end
                                
                                % compute the no of constraints satisfied
                                % by this split.                                
                                nsatis_temp = sum(clusters_m_temp(CL_m(:,1)) ~= clusters_m_temp(CL_m(:,2)));
                                %Check if we're better
                                if cut_m1_temp + cut_m2_temp < cut_m1 + cut_m2                                    
                                    [cut_m1, cut_m2, clusters_m, nsatis]=deal(cut_m1_temp, cut_m2_temp, clusters_m_temp, nsatis_temp);
                                    %assert(length(clusters_m)==length(ix_m));
                                end
                                %assert(logical(exist('clusters_m','var')));
                            end
                        end
                    catch err
                        if verbosity >= 3
                            display(err);
                        end
                        cut_m1 = inf;
                        cut_m2 = inf;
                        nsatis = 0;
                    end
                end
                %if(true)
                if(cut_m1+cut_m2>0)
                    if verbosity >= 2
                        display('The isolated components have non-zero cut to the rest of the graph');
                    end
                    for m1=1:length(sizes)
                        ix_comp=find(comp==m1);
                        % if the size of the current connected component is larger than 1, try to partition it
                        if (length(ix_comp)>1)
                            Wm_comp=sparse(Wm(ix_comp,ix_comp));
                            Wm_comp2=Wm_comp;
                            scaleOfW = 1;
                            if(2*max(sum(Wm_comp.^2))<eps)
                                Wm_comp2=Wm_comp/max(max(Wm_comp));
                                scaleOfW = max(max(Wm_comp));
                            end
                            
                            
                            vertex_weights_m_comp = vertex_weights_m(ix_comp);
                            Q_m_comp = Q_m(ix_comp, ix_comp);
                            [cix, cjx] = find(triu(Q_m_comp));
                            CL_m_comp = [cix cjx];
                            if(~connected && verbosity>=3), fprintf('%2sComputing partitioning of connected component %d of %d of subgraph %d.\n',' ', m1,length(sizes), m); end
                            if (~connected)
                                ix_rest=find(comp~=m1); % all other components in the current cluster
                                cut_rest=sum(sum(W(ix_m(ix_rest),setdiff(1:n,ix_m(ix_rest)))));
                                size_rest=sum(vertex_weights(ix_m(ix_rest)));
                            else
                                cut_rest=0;
                                size_rest=0;
                            end
                            
                            result_ix = result_ix_conn_comp(min(ix_m(ix_comp)), length(ix_comp));
                            if verbosity >= 2
                                display(['Trying to retrieve the previously computed result from the index: (', num2str(min(ix_m(ix_comp))), ', ', num2str(length(ix_comp)), ')', ' and the location: ', num2str(result_ix)]);
                            end
                            if result_ix~= 0
                                vmin_m_temp = result_conn_comp_vmin{result_ix};
                                [clusters_m_temp, cut_m1_temp, cut_m2_temp] = opt_thresh_cnstr_functional_subgraph(vmin_m_temp, Wm_comp, deg(ix_m(ix_comp)), vertex_weights_m_comp', Q_m_comp,0, 1, ix_comp, ix_m, cut_rest, size_rest, size_m);
                                if verbosity >= 2
                                    display(['Retrieved the previously computed result from the index: (', num2str(min(ix_m(ix_comp))), ', ', num2str(length(ix_comp)), ')', ' and the location: ', num2str(result_ix)]);                                
                                end
                            else
                                
                                
                                %cnstr_prob_start = tic;
                                %min_ncl = 0;
                                knew = size(Wm_comp2,1)*k/n;
                                knew = max(2, floor(knew));
                                %min_ncl = 0.5*(knew - 2)*size(CL,1)/(knew-1);
                                %min_ncl = floor(min_ncl);
                                min_ncl = (knew - 2)*size(CL_m,1)/knew;
                                min_ncl = floor(min_ncl);
                                k1 = floor(knew/2);
                                k2 = knew-k1;
                                max_ncl = size(CL_m,1)* (k1*(k1-1)+k2*(k2-1))/(knew *(knew-1));
                                max_ncl = floor(max_ncl);
                                max_ncl = min_ncl;
                                %[vmin,cut, clusters, vmins, cuts_all, objectives, viols]= ...
                                %if verbosity >=1, fprintf('Minimizing the Functional using %d starting points\n', nRuns+start_flags); end
                                if verbosity >=1
                                    if ~isempty(CL), fprintf('%2sMinimizing the Functional subject to Cannot link constraints using %d starting points\n', ' ', nRuns+start_flags);
                                    else fprintf('%4sMinimizing the unconstrained Functional using %d starting points\n', ' ', nRuns+start_flags);
                                    end
                                end
                                cnstr_prob_start = tic;
                                [vmin_comp, cut_m1_temp, cut_m2_temp, clusters_m_temp] = ...
                                    hierarchical_solve_cnstr_functional_incremental_subgraph(Wm_comp2,CL_m_comp, deg(ix_m(ix_comp)),vertex_weights_m_comp, ...
                                    start_flags, nRuns, min_ncl, max_ncl, MAX_RUNS, FISTA,MAX_ITERS,verbosity, ix_comp', ix_m', cut_rest, size_rest, size_m, scaleOfW);
                                %                                     hierarchical_solve_cnstr_functional_incremental_subgraph(Wm_comp2,CL_comp,vertex_weights_comp,start_flags, nRuns, min_ncl, FISTA,verbosity>=3);
                                %                                 start_comp,Wm_comp,normalized,deg,criterion_threshold,index_comp,index_m,cut_rest,size_rest,size_m);
                                
                                %toc(cnstr_prob_start);                                    
                                
                                if(verbosity>=2)
                                    displayCurrentObjective(sum(cutParts)-cutParts(m)+cut_m1_temp+cut_m2_temp,-1,true);
                                    display(['Time taken for one split: ', num2str(toc(cnstr_prob_start))]);                                
                                end
                                % Save this result if it is one of the several connected components of the subgraph
                                if length(sizes)>1
                                    connectedComponentCounter = connectedComponentCounter+1;
                                    result_conn_comp_vmin{connectedComponentCounter} = vmin_comp;
                                    result_ix_conn_comp(min(ix_m(ix_comp)), length(ix_comp)) = connectedComponentCounter;
                                    if verbosity >= 2
                                        display(['Saved the computed result to the index: (', num2str(min(ix_m(ix_comp))), ', ', num2str(length(ix_comp)), ')', ' and the location: ', num2str(connectedComponentCounter)]);
                                    end
                                end
                            end
                            
                            nsatis_temp = sum(clusters_m_temp(CL_m(:,1)) ~= clusters_m_temp(CL_m(:,2)));
                            % Check if we're better
                            if cut_m1_temp + cut_m2_temp < cut_m1 + cut_m2
                                [cut_m1, cut_m2, clusters_m, nsatis]=deal(cut_m1_temp, cut_m2_temp, clusters_m_temp, nsatis_temp);
                                %assert(length(clusters_m)==length(ix_m));
                            end
                        end
                    end
                end
                % store current best partition
                sub_clusters{m}=clusters_m;
                sub_cuts(m,1) = cut_m1;
                sub_cuts(m,2) = cut_m2;
                sub_nsatis(m) = nsatis;
                nCL_m = size(CL_m,1);
                sub_nCL(m) = nCL_m;
                
            end
            
            % print out best cut possible by partitioning of current subgraph
            cut = sum(cutParts)-cutParts(m)+cut_m1+cut_m2;
            if(verbosity>=2)
                fprintf('\n%4sBest result achievable by partitioning of subgraph %d:\n',' ', m);
                displayCurrentObjective(cut,-1,true);
                fprintf('\n');
            end
            
            nviols = ncl - nsatis_tot - nsatis;
            % check if partitoning of the current subgraph gives better cut
            %if ((nCL_m>0) || best_cut==inf) && cut+nviols/ncl<best_cut+ (ncl - nsatis_tot - best_satis)/ncl
            %if cut<best_cut
            %if nviols/ncl<(ncl - nsatis_tot - best_satis)/ncl
            
            flag = false;
            knew = size(clusters_m,1)*k/n; knew = max(2, floor(knew));
            exp_cls = knew*(knew-1)*size(CL,1)/(k*(k-1));
            min_satis = 2*exp_cls/knew;
            min_satis = floor(min_satis);
            
            %if best_cut == inf
            %   flag = true;
            
            %                    else
            
            if nsatis > current_best_nsatis % this is not used currently
                current_best_nsatis = nsatis;
                current_best_m = m;
            end
            
            % We want a certain minimum number of constraints to be
            % satisfied in each round and an expected total number until
            % the current round.
            if nsatis >= min_satis
                
                if ~isempty(CL) && l==k % In the final iteration, we want to satisfy as many constraints as possible.
                    %nsatis_tot >= l*(l-1)*size(CL,1)/(k*(k-1)) || l==k
                    if nsatis > best_satis % we choose the split that satisfies many constraints.
                        flag = true;
                    else if nsatis == best_satis && cut < best_cut % among those splits that satisfy the same number of constraints we choose the one that has better cut.
                            flag = true;
                        else
                            flag = false;
                        end
                    end
                    
                    %nsatis >=best_satis || nsatis + nsatis_tot >= l*(l-1)*size(CL,1)/(k*(k-1)) % it should be (2k-l)*(l-1) instead of l(l-1)!
                    %reasoning: after the current split only the
                    %points in the remaining (n-k+1) classes
                    %violate the constraints: (n-k+1)*(n-k)/2. The
                    %total satis is obtained by subtracting above
                    %from k*(k-1)/2
                % best satis is the maximum number of constraints satisfied
                % in this round among those splits that satisfied min_satis
                % required.
                % nsatis_tot is the total number of constraints satisfied
                % until the last round.
                % We will compare the cut value of those splits that have
                % hightest number of constraints satisfied or that have
                % satisfied in total an expected number until this round.
                elseif nsatis >= best_satis || nsatis + nsatis_tot >= (2*k-l)*(l-1)*size(CL,1)/(k*(k-1)) % it should be (2k-l)*(l-1) instead of l(l-1)!
                    
                    if cut<best_cut
                        flag=true;
                    end
                end
            end
            
            if flag
                [best_cut,bestCheeger,best_cut_m1,best_cut_m2,best_satis,best_m]= deal(cut,-1,cut_m1,cut_m2,nsatis,m);
                clusters_new=clusters;
                clusters_new(ix_m)=(l-m)*clusters_m+clusters_new(ix_m);
                
                %assert(best_cut>=0);% && bestCheeger>=0);
            end
            
            %                 if best_cut == inf || nsatis >= best_satis || nsatis + nsatis_tot >= l*(l-1)*size(CL,1)/(k*(k-1));
            %                     if cut < best_cut
            %                         [best_cut,bestCheeger,best_cut_m1,best_cut_m2,best_satis,best_m]= deal(cut,-1,cut_m1,cut_m2,nsatis,m);
            %                         clusters_new=clusters;
            %                         clusters_new(ix_m)=(l-m)*clusters_m+clusters_new(ix_m);
            %
            %                         assert(best_cut>=0);% && bestCheeger>=0);
            %                     end
            %                 end
            
            % if we have already found a partition with cut 0, we don't
            % need to consider the other subgraphs
            if best_cut==0
                break;
            end
        end
        
        %knew = floor(knew/2);
        
        if(best_cut==inf) % i.e. none of them were selected as none satisfied min expected (or it is the final iteration)... so choose the one that satisfies many
            
            max_satis = max(sub_nsatis(1:l-1));
            ixs = find(sub_nsatis(1:l-1) == max_satis);
            m = ixs(1);
            cut_best = sum(cutParts)-cutParts(m)+sub_cuts(m,1)+sub_cuts(m,2);
            
            for jj=2:length(ixs)
                cut_cur = sum(cutParts)-cutParts(ixs(jj))+sub_cuts(ixs(jj),1)+sub_cuts(ixs(jj),2);
                if cut_cur < cut_best
                    
                    m = ixs(jj);
                    cut_best = cut_cur;
                end
            end
            
            %m = current_best_m;
            ix_m = find(clusters==m);
            clusters_m = sub_clusters{m};
            cut_m1 = sub_cuts(m,1);
            cut_m2 = sub_cuts(m,2);
            nsatis = sub_nsatis(m);
            nCL_m = sub_nCL(m);
            
            cut = sum(cutParts)-cutParts(m)+cut_m1+cut_m2;
            [best_cut,bestCheeger,best_cut_m1,best_cut_m2,best_satis,best_m]= deal(cut,-1,cut_m1,cut_m2,nsatis,m);
            clusters_new=clusters;
            if ~isempty(clusters_m)
                clusters_new(ix_m)=(l-m)*clusters_m+clusters_new(ix_m);
            else
                clusters_new(ix_m)=clusters_new(ix_m);
            end
            
            %assert(best_cut>=0);
        end
        
        % Update
        nsatis_tot = nsatis_tot + best_satis;
        clusters=clusters_new;
        cuts_intermediate(1,l-1)=best_cut;
        clusters_intermediate(:,l-1)=clusters;
        
        cutParts(best_m)=best_cut_m1;
        cutParts(l)=best_cut_m2;
        
        % Check that we have the right number of clusters
        %assert(length(unique(clusters))==l);
        
        % Reset subcutparts and subclusters;
        sub_cuts(best_m,:)=0;
        sub_clusters{best_m}= [];
        sub_cuts(l,:)=0;
        sub_clusters{l}= [];
        
        % Print out current objective
        if(verbosity>=1)
            fprintf('%2sDecided to partition subgraph %d. Finished Clustering into %d parts.\n',' ', best_m,l);
            displayCurrentObjective(best_cut,-1,true);
            fprintf('\n');
        end
        
    end
    cut = cuts_intermediate(end);
    viols = 0;
    if ~isempty(CL)
        viols = sum(clusters(CL(:,1)) == clusters(CL(:,2)));
    end
    
    
    clusters = clusters(map);
    
    for i=1:k-1
        cl = clusters_intermediate(:,i);
        clusters_intermediate_orig(:,i) = cl(map);
    end
    
end

end


function displayCurrentObjective(cut_temp,cheeger_temp,normalized)

if (normalized)
    fprintf('...Normalized Cut: %.8g   Normalized Cheeger Cut: %.8g\n',cut_temp,cheeger_temp);
else
    fprintf('...Ratio Cut: %.8g   Ratio Cheeger Cut: %.8g\n',cut_temp,cheeger_temp);
end

end
