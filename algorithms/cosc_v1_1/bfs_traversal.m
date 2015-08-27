function [labels] = bfs_traversal(W)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    n = size(W,1);
    visited = false(n,1);
    labels = zeros(n,1);
    
    queue = 1;
    labels(1) = 1;
    visited(1) = true;
    while ~isempty(queue)
    
        v = queue(1);
        queue = queue(2:end);
        
	    neigh = find(W(v,:));
        unlabeled_neigh = neigh( ~visited(neigh) );
        labels( unlabeled_neigh  ) = ~labels(v);
        
%         if sum( labels(setdiff(neigh, unlabeled_neigh)) == labels(v) ) > 0 
%             display('inconsistency');
%         end
        
        queue = [queue unlabeled_neigh];
        visited( unlabeled_neigh ) = 1;
        
        
    end
    
    assert( sum(visited) == n );
    
end