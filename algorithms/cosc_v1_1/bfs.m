function comp = bfs(W, u, n)
% Usage comp = bfs(W, u, n)
% u is starting vertex, n is the number of vertices.
% 
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    visited = sparse(n,1);
%     visited(u) = 1;
    
    unexplored = u;
    visited( unexplored ) = 1;  % one step after visited is set to 1, the corresponding vertices are explored.
    while ~isempty(unexplored)
    
        [neigh, jj] = find(W(:,unexplored));   
        unexplored = neigh( ~visited(neigh) );
        visited( unexplored ) = 1;

    end
    
    comp = find(visited);
    
    % Sanity check
    [ix, jx] = find(W(:,comp));
    ix = [ix; u]; % add the starting vertex!
    assert( isempty(setdiff(comp, ix)) );  % We have not missed any vertex that is reachable from current component.

end
