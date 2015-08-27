function [clusters,cuts,cheegers, vmin, W, vertex_weights, Y, cnstr1, cnstr2] = ...
    hierarchical_spec_clustering_lcnstrs_merging(W, vertex_weights, CL, ML, normalized, verbosity, Y)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    cnstr1 = [];
    cnstr2 = [];
    
    if ~isempty(CL)
        cnstr1 = CL(:,1);
        cnstr2 = CL(:,2);
    end
    
    [W, vertex_weights, map, Y] = process_mls( W, vertex_weights, ML, Y );
    cnstr1 = map(cnstr1);
    cnstr2 = map(cnstr2);   
    
    W = triu(W);
    W = W + W';
    
    [clusters,cuts,cheegers, vmin, W, vertex_weights, Y, cnstr1, cnstr2] = ...
        spec_clustering_lcnstrs(W, vertex_weights, [cnstr1 cnstr2], [], normalized, verbosity, Y);

    clusters = clusters(map_dml);
    clusters = clusters(map);

end