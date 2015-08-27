function [color, comp] = two_coloring(W, u, n)
% Usage color = two_coloring(W, u, n)
% u is starting vertex, n is the number of vertices.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
    visited = sparse(n,1);
    color = false(n,1);
    
    unexplored = u;
    visited( unexplored ) = 1;  % one step after visited is set to 1, the corresponding vertices are explored.
    color(unexplored) = true;
    while ~isempty(unexplored)
        
        current_color = color(unexplored(1));
        [neigh, jj] = find(W(:,unexplored));   
        unexplored = neigh( ~visited(neigh) );
        color(unexplored) = ~current_color;
        visited( unexplored ) = 1;

    end

    % Sanity check
    comp = find(visited);
    assert( sum(visited) == length(comp) );  % every vertex colored
    [ix, jx] = find(W(comp, comp));     
 %   assert( sum(color(comp(ix)) == color(comp(jx))) == 0 );   % and colors respect given edges
    
%    display(sum(color(comp(ix)) == color(comp(jx))) == 0 );   % and colors respect given edges
end
