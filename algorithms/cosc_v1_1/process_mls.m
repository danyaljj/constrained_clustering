function [Wp, vertex_weights_p, map,prev_clustersp] = process_mls( W, vertex_weights, ML, prev_clusters)
% Usage: [Wp, vertex_weights_p, map] = process_mls( W, vertex_weights, ML)
%
%   Output:
%           Wp = Weight matrix of the sparsified graph
%           vertex_weights_p = vertex weights of the sparsified graph
%           map = map to the new indices in the sparsified graph
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    


    Wp = W;
    vertex_weights_p = vertex_weights;
    map = 1:size(W,1);
    map = map';
    prev_clustersp = prev_clusters;

    if ~isempty(ML)
        
        ml_vtxs = ML(:);
        ml_vtxs = unique(ml_vtxs);        
        
        cur_ix = 1;
        u = ml_vtxs(cur_ix);
        n = size(W,1);

        ixs = cell(n,1);
        ilarge_comp = 0;
        
        visited = sparse(size(W,1),1);
        %unvisited = ones(size(W,1),1);
        
        % Construct constraint graph
        Q = construct_cnstr_graph( ML, n );

        % start breadth first search on Q from u.
        comp = bfs(Q, u, n);
%         comp1 = graphtraverse(Q, u, 'METHOD', 'BFS', 'DIRECTED', false); %comp = comp1';
%         assert( sum( sort(comp1) - comp ) == 0 );
        visited(comp) = 1;
%        unvisited(comp)=0;
        if length(comp) > 1
            ilarge_comp = ilarge_comp + 1;
            ixs{ilarge_comp} = comp;
        end

        
        %while sum(visited) < n
        while cur_ix < length(ml_vtxs)

            %u = find(visited==0, 1);
            %u = find(unvisited, 1);
            cur_ix = cur_ix+1;
            u = ml_vtxs(cur_ix);

            while visited(u)
               cur_ix = cur_ix+1;
               if cur_ix > length(ml_vtxs)
                   break;
               end
               u = ml_vtxs(cur_ix);                
            end
            
            if cur_ix > length(ml_vtxs)
                   break;
            end
            comp = bfs(Q, u, n);
%             comp1 = graphtraverse(Q, u, 'METHOD', 'BFS', 'DIRECTED', false);%comp = comp1';
%             assert( sum( sort(comp1') - comp ) == 0 );
            visited(comp) = 1;
            %unvisited(comp) = 0;
            if length(comp) > 1
                ilarge_comp = ilarge_comp + 1;
                ixs{ilarge_comp} = comp;
            end                

        end

        [Wp, vertex_weights_p, map, prev_clustersp] = merge( W, vertex_weights, ixs(1:ilarge_comp), prev_clusters);
    end

end