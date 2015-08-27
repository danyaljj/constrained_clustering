%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%    optimal thresholding of our constrained functional on a subgraph
%%    (possibly several disconnected components)
%%    numerator is half the cut + half the 
%%
%%    _m refers to the whole subgraph
%%    _comp refers to the component on which constrained 1-spectral
%%    clustering mehtod is run. 
%%    _rest refers to the remaining disconnected components. 
%%
%%    the input argument index_m is not requied/used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [allClustersInClusterM, cutPart1,cutPart2, cnstrpart] = ...
    opt_thresh_cnstr_functional_subgraph(vmin_comp, W, deg_comp, vertex_weights_comp,Q, gamma, criterion_threshold, ...
    index_comp,index_m,cut_rest,size_rest,size_m)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%          
        
        if (size(deg_comp,1)>1) 
            deg_comp=deg_comp';
        end
        
        %% Qvol is the total number of cannot link constraints on the
        %% subgraph
        Qvol = sum(sum(Q));
        Qvol = Qvol/2;
        
        
%         gamma = gamma/2;
%         deg=full(sum(W));

        cnstr_deg = full(sum(Q));
%         W_comp = W(index_comp, index_comp);
%         Q_comp = Q(index_comp, index_comp);

        W_comp = W;
        Q_comp = Q;
        
        [vminM_sorted, index]=sort(vmin_comp);
        [vminU,indexU]=unique(vminM_sorted);
%         vminU = vminM_sorted; indexU = 1:size(W,1);
        
        
        W_sorted=W_comp(index,index);
        Q_sorted=Q_comp(index,index);
        vertex_weights_comp = vertex_weights_comp(index);
        balance_thresholds = cumsum(vertex_weights_comp);

        % calculate cuts
        %deg_comp=deg(index_m(index_comp));
        %cnstr_deg_comp=cnstr_deg(index_m(index_comp));
        cnstr_deg_comp=cnstr_deg;
        volumes_threshold=cumsum(deg_comp(index));
        cnstr_volumes_threshold=cumsum(cnstr_deg_comp(index));
        triup=triu(W_sorted);
        Q_triup=triu(Q_sorted);
        tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triup)));
        cnstr_tempcuts_threshold=cnstr_volumes_threshold - 2*cumsum(full(sum(Q_triup)));
        tempcuts_threshold2=(volumes_threshold(end)-volumes_threshold) - (sum(sum(W_sorted))-2*cumsum(full(sum(triup,2)))');            
        cnstr_tempcuts_threshold2=(cnstr_volumes_threshold(end)-cnstr_volumes_threshold) - (sum(sum(Q_sorted))-2*cumsum(full(sum(Q_triup,2)))');            

        tempcuts_threshold=tempcuts_threshold(indexU);
        cnstr_tempcuts_threshold=cnstr_tempcuts_threshold(indexU);
        tempcuts_threshold2=tempcuts_threshold2(indexU);
        cnstr_tempcuts_threshold2=cnstr_tempcuts_threshold2(indexU);

        %volumes_threshold=volumes_threshold(indexU); %balance_thresholds will be used instead of volumes_thresholds
        balance_thresholds = balance_thresholds(indexU);
        
        
        % divide by size/volume
        %if(normalized)
        % Option 1: merge the remaining components and C_t.
        % In this case, we need to count the constraints voilated
        % between the components C_t^prime and the merged part.
        cutparts1_threshold=(tempcuts_threshold(1:end-1)+cut_rest)./(balance_thresholds(1:end-1)+size_rest);
        cnstrparts1_threshold=gamma*(Qvol - cnstr_tempcuts_threshold(1:end-1))./(balance_thresholds(1:end-1));
        cutparts2_threshold=tempcuts_threshold2(1:end-1)./(balance_thresholds(end)-balance_thresholds(1:end-1));
        cnstrparts2_threshold=gamma*(Qvol - cnstr_tempcuts_threshold2(1:end-1))./(balance_thresholds(end)-balance_thresholds(1:end-1));

        % Option 2: merge the remaining components and C_t^prime
        % In this case, we need to count the constraints voilated
        % between the components C_t and the merged part.
        cutparts1b_threshold=tempcuts_threshold(1:end-1)./balance_thresholds(1:end-1);
        cnstrparts1b_threshold=gamma*(Qvol-cnstr_tempcuts_threshold(1:end-1))./balance_thresholds(1:end-1);
        cutparts2b_threshold=(tempcuts_threshold2(1:end-1)+cut_rest)./((balance_thresholds(end)-balance_thresholds(1:end-1))+size_rest);
        cnstrparts2b_threshold=gamma*(Qvol-cnstr_tempcuts_threshold2(1:end-1))./((balance_thresholds(end)-balance_thresholds(1:end-1)));
            
        
        %calculate total cuts
        if(criterion_threshold==1)

            cuts_threshold=cutparts1_threshold+cutparts2_threshold + cnstrparts1_threshold+cnstrparts2_threshold;
            [cut1,threshold_index]=min(cuts_threshold);
            
            
            cuts_threshold_b=cutparts1b_threshold+cutparts2b_threshold + cnstrparts1b_threshold+cnstrparts2b_threshold;
            [cut1b,threshold_index_b]=min(cuts_threshold_b);
                        
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
            cnstrpart = cnstrparts1_threshold(threshold_index)+cnstrparts2_threshold(threshold_index); 
            
            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index));
        
            allClustersInClusterM= zeros(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        else
            cutPart1=cutparts1b_threshold(threshold_index_b);
            cutPart2=cutparts2b_threshold(threshold_index_b);
            cnstrpart = cnstrparts1_threshold(threshold_index_b)+cnstrparts2_threshold(threshold_index_b);
            allClustersInClusterM_comp= (vmin_comp>vminU(threshold_index_b));
        
            allClustersInClusterM= ones(size_m,1);
            allClustersInClusterM(index_comp)=allClustersInClusterM_comp;
        end
        
%    cuts_threshold

end

