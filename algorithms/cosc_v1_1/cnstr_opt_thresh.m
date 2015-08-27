function [allClustersInClusterM, cutPart1,cutPart2, cnstrPart, cnstr_Sat_m, hbcut, hCutPart1, hCutPart2] = ...
    cnstr_opt_thresh(vmin_comp,W_comp,W,normalized, vertex_weights_comp, ...
    criterion_threshold,index_comp,index_m,cut_rest,size_rest,size_m, Qm, Qm_comp, lambda, totCnstrs, hbcut, Y, flag)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
      
        deg=full(sum(W));
        
        [vminM_sorted, index]=sort(vmin_comp);
        [vminU,indexU]=unique(vminM_sorted);
%          vminU = vminM_sorted;
%          indexU = 1:size(W,1);
        
        W_sorted=W_comp(index,index);
        Qm_comp_sorted=Qm_comp(index,index);
        vertex_weights_comp = vertex_weights_comp(index);
        balance_thresholds = cumsum(vertex_weights_comp);

        % calculate cuts
        deg_comp=deg(index_m(index_comp));
        volumes_threshold=cumsum(deg_comp(index));
        triup=triu(W_sorted);
        tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triup)));
        tempcuts_threshold2=(volumes_threshold(end)-volumes_threshold) - (sum(sum(W_sorted))-2*cumsum(full(sum(triup,2)))');            

        %calculate fraction voilated for each threshold:
        % cut computed on the matrix Q, gives us the number of satisfied
        % constraints

        % Consider the constraint matrix (Q) on the whole subgraph.
        % Sort the current component according to vmin_comp
        % Append the remaining connected components to the end 
        Q_index_m = [index_comp(index); setdiff((1:length(index_m))', index_comp)];
        %Qdm = full(sum(Qm));
        %Q_volumes_threshold=cumsum(Qdm(Q_index_m));
        Qm_comp_rectangular = Qm( index_comp(index), Q_index_m );
        Qdm = full(sum(Qm_comp_rectangular,2));
        Q_volumes_threshold=cumsum(Qdm)';
        
        % assert Qm(Q_index_m, Q_index_m) and Qm_comp_sorted are same!
        Q_triup=triu(Qm_comp_sorted);
        Q_triup1=triu( Qm(index_comp(index),index_comp(index)) );
        assert( sum(sum(Q_triup - Q_triup1)) == 0 );
        
        % Now there are two options: (1) merge the remaining components to
        % the component C_t^prime.
        Q_tempcuts_threshold1_temp=Q_volumes_threshold - 2*cumsum(full(sum(Q_triup)));

        %(2) merge the remaining components to the component C_t
        Q_tempcuts_threshold2_temp=(Q_volumes_threshold(end)-Q_volumes_threshold) - (sum(sum(Qm_comp_sorted))-2*cumsum(full(sum(Q_triup,2)))');
        %The above gives the number of satisfied constraints... now obtain fraction
        %voilated..
        if ~isempty(find(Qm,1))
            Q_tempcuts_threshold1=(0.5*sum(sum(Qm)) - Q_tempcuts_threshold1_temp)./totCnstrs;
            Q_tempcuts_threshold2=(0.5*sum(sum(Qm)) - Q_tempcuts_threshold2_temp)./totCnstrs;
       else
            Q_tempcuts_threshold1=Q_tempcuts_threshold1_temp;
            Q_tempcuts_threshold2=Q_tempcuts_threshold2_temp;
       end
        
        tempcuts_threshold=tempcuts_threshold(indexU);
        tempcuts_threshold2=tempcuts_threshold2(indexU);
        volumes_threshold=volumes_threshold(indexU);
        balance_thresholds = balance_thresholds(indexU);
        
        Q_tempcuts_threshold1=Q_tempcuts_threshold1(indexU);
        Q_tempcuts_threshold2=Q_tempcuts_threshold2(indexU);
        
        
        % divide by size/volume
        if(normalized)
            % Option 1: merge the remaining components and C_t.
            % In this case, we need to count the constraints voilated
            % between the components C_t^prime and the merged part.
            cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(balance_thresholds(1:end-1)+size_rest);
            cnstrparts1_threshold=Q_tempcuts_threshold2(1:end-1);
            
            cutparts2_threshold=tempcuts_threshold2(1:end-1)./(balance_thresholds(end)-balance_thresholds(1:end-1));
            %cnstrparts2_threshold=lambda*Q_tempcuts_threshold2(1:end-1);

            % Option 2: merge the remaining components adn C_t^prime
            % In this case, we need to count the constraints voilated
            % between the components C_t and the merged part.
            cutparts1b_threshold=tempcuts_threshold(1:end-1)./balance_thresholds(1:end-1);
            cnstrparts1b_threshold=Q_tempcuts_threshold1(1:end-1);
            
            cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((balance_thresholds(end)-balance_thresholds(1:end-1))+size_rest);
            %cutparts2b_threshold=cutparts2b_threshold + lambda*Q_tempcuts_threshold1(1:end-1);
        else
            sizes_threshold=cumsum(ones(1,size(vmin_comp,1)-1));
            sizes_threshold=sizes_threshold(indexU(1:end-1));
            cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(sizes_threshold+size_rest);
            cnstrparts1_threshold=Q_tempcuts_threshold2(1:end-1);
            
            cutparts2_threshold=tempcuts_threshold2(1:end-1)./(size(vmin_comp,1)-sizes_threshold);
            %cutparts2_threshold=cutparts2_threshold + lambda*Q_tempcuts_threshold2(1:end-1);
            
            cutparts1b_threshold=tempcuts_threshold(1:end-1)./sizes_threshold;
            cnstrparts1b_threshold=Q_tempcuts_threshold1(1:end-1);
            cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((size(vmin_comp,1)-sizes_threshold)+size_rest);
            %cutparts2b_threshold=cutparts2b_threshold + lambda*Q_tempcuts_threshold1(1:end-1);
        end

        
        
        %calculate total cuts
        if(criterion_threshold==1)
            [hbcut, h_index] = min(cutparts1_threshold+cutparts2_threshold);
            if lambda == inf
                min_voilated = min(cnstrparts1_threshold);
                min_indices = find(cnstrparts1_threshold == min_voilated);
                cuts_threshold=(cutparts1_threshold+cutparts2_threshold);
                [cut1, threshold_index] = min(cuts_threshold(cnstrparts1_threshold == min_voilated));
                threshold_index = min_indices(threshold_index);
            else
                cuts_threshold=(cutparts1_threshold+cutparts2_threshold)/hbcut +lambda*cnstrparts1_threshold;
                [cut1,threshold_index]=min(cuts_threshold);
            end
            
            
            [hbcutb, h_indexb] = min(cutparts1b_threshold+cutparts2b_threshold);
            if lambda == inf
                min_voilated = min(cnstrparts1b_threshold);
                min_indices = find(cnstrparts1b_threshold == min_voilated);
                cuts_thresholdb=(cutparts1b_threshold+cutparts2b_threshold);
                [cut1b, threshold_index_b] = min(cuts_thresholdb(cnstrparts1b_threshold == min_voilated));
                threshold_index_b = min_indices(threshold_index_b);
            else
                cuts_threshold_b=(cutparts1b_threshold+cutparts2b_threshold)/hbcutb +lambda*cnstrparts1b_threshold;
                [cut1b,threshold_index_b]=min(cuts_threshold_b);
            end
% 
%             [hbcut, h_index] = min(cutparts1_threshold+cutparts2_threshold);
%             cuts_threshold=(cutparts1_threshold+cutparts2_threshold)/hbcut +lambda*cnstrparts1_threshold;
%             [cut1,threshold_index]=min(cuts_threshold);
%             
%             [hbcutb, h_indexb] = min(cutparts1b_threshold+cutparts2b_threshold);
%             cuts_threshold_b=(cutparts1b_threshold+cutparts2b_threshold)/hbcutb +lambda*cnstrparts1b_threshold;
%             [cut1b,threshold_index_b]=min(cuts_threshold_b);
            
            comp_case=1;
            if (cut1b<cut1) 
                comp_case=2;
            end
        else
            cheegers_threshold=max(cutparts1_threshold,cutparts2_threshold);
            [cheeger1,threshold_index]=min(cheegers_threshold);
            
            cheegers_threshold_b=max(cutparts1_threshold,cutparts2_threshold);
            [cheeger1b,threshold_index_b]=min(cheegers_threshold_b);
            
            comp_case=1;
            if (cheeger1b<cheeger1) 
                comp_case=2;
            end
        end

        if(comp_case==1)
            cutPart1=cutparts1_threshold(threshold_index);
            cutPart2=cutparts2_threshold(threshold_index);
            
            hCutPart1 = cutparts1_threshold(h_index);
            hCutPart2 = cutparts2_threshold(h_index);
            
            cnstr_Sat_m = 0.5*sum(sum(Qm)) - cnstrparts1_threshold(threshold_index) * totCnstrs;
            cnstrPart = lambda*cnstrparts1_threshold(threshold_index);

            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index));
        
            allClustersInClusterM= zeros(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        else
            cutPart1=cutparts1b_threshold(threshold_index_b);
            cutPart2=cutparts2b_threshold(threshold_index_b);
            
            hCutPart1 = cutparts1_threshold(h_indexb);
            hCutPart2 = cutparts2_threshold(h_indexb);
            
            cnstr_Sat_m = 0.5*sum(sum(Qm)) - cnstrparts1b_threshold(threshold_index) * totCnstrs;
            cnstrPart = lambda*cnstrparts1b_threshold(threshold_index);

            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index_b));
        
            allClustersInClusterM= ones(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        end
        
    %cuts_threshold
    
    if flag && min_voilated > 0
        allClustersInClusterM = inf;
    end
    
%     disp(thresho

end

