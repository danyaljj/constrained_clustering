function [Wp, bp, map, prev_clustersp] = merge( W, b, ixs, prev_clusters )
% Usage: [Wp, b, map] = sparsify( W, b, ix )
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    n = size(W,1);
    k = 0;
    ixs_all = [];
    for i=1:length(ixs)
        k = k+length(ixs{i}) - 1;
        ixs_all = [ixs_all; ixs{i}];
    end
    
    ixsc = setdiff(1:n, ixs_all);
    np = n - k;
    n_ixsc = length(ixsc);
    
    Wp = sparse(np, np);

    bp = zeros(np, 1);
    bp(1:n_ixsc) = b(ixsc);
    
    prev_clustersp = cell(size(prev_clusters));
    for i=1:length(prev_clustersp)
        prev_clustersp{i} = inf*ones(np,1);
        prev_clustersp{i}(1:n_ixsc) = prev_clusters{i}(ixsc);
    end        

    map = zeros(n,1);
    map(ixsc) = 1:n_ixsc;
    for i=1:length(ixs)
        map(ixs{i}) = n_ixsc + i;
    end      

    all_indices = 1:n;
    all_indices_p = map(all_indices);
    %temp_W = sparse(np, length(ixs));
    for i=1:length(ixs)
        
        subset = ixs{i};
%        bp(n_ixsc+i) = sum(b(subset));
        tot_weight = 0;
        for index=1:length(subset)
            tot_weight = tot_weight + b(subset(index));
        end
        bp(n_ixsc+i) = tot_weight;        
        
        for jj=1:length(prev_clustersp)
            % Choose the major label of the vertices being merged 
            % This makes no sense for a random start!
%             labels = unique(prev_clusters{jj});
%             nix1 = sum(prev_clusters{jj}(subset)==labels(1));
%             major_label = nix1;
% 
%             if length(labels) > 1
%                 nix2 = sum(prev_clusters{jj}(subset)==labels(2));
% 
%                 if nix2 > nix1
%                     major_label = nix2;
%                 end
%             end
                
            prev_clustersp{jj}(n_ixsc+i) = prev_clusters{jj}(subset(1));
        end
                
            
        Wp(:, n_ixsc+i) = sparse( all_indices_p,1, sum( W(:, subset), 2 ), np, 1);
        %temp_W(:,i) = sparse( all_indices_p,1, sum( W(:, subset), 2 ), np, 1);
        
    end    
    
    %Wp(:, n_ixsc+1:n_ixsc+length(ixs)) = temp_W;
    Wp = [W(ixsc, ixsc) Wp(1:n_ixsc, n_ixsc+1:end); Wp(:, n_ixsc+1:end)'];
   
    %set the diagonal entries to zero
    Wp(1:np+1:np*np) = 0;
end