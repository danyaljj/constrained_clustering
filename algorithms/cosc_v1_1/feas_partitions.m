function fpartitions = feas_partitions( W, CL, Y, C )
%% the kth feasible solution, fsolns(:, k) has three values: 0, 1, -10. 
%% where -10 is for unconstrained vertices.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    n = size(W,1);
    fpartitions = -10*ones(n,1);
    
    if ~isempty(CL)
        
        u = 1;
%         u = CL(1,1);
        visited = sparse(size(W,1),1);
        
        % Construct constraint graph
        Q = construct_cnstr_graph( CL, n );

        % start breadth first search on Q from u.
        comp = bfs(Q, u, n);
        visited(comp) = 1;
        if length(comp) > 1
            Qcomp = Q(comp, comp);
            color = two_coloring(Qcomp, 1, size(Qcomp,1));
            fpartitions(comp) = color;
        end

        
        while sum(visited) < n

            u = find(visited==0, 1);
            comp = bfs(Q, u, n);
            visited(comp) = 1;
            if length(comp) > 1
                Qcomp = Q(comp, comp);
                color = two_coloring(Qcomp, 1, size(Qcomp,1));
                fpartitions(comp) = color;
            end                

        end
    end

   assert( sum(fpartitions(CL(:,1)) == fpartitions(CL(:,2))) == 0 );
   
end