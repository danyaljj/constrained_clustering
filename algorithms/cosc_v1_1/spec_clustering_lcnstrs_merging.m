function [clusters,cuts,cheegers, vmin, W, vertex_weights, Y, cnstr1, cnstr2] = ...
    spec_clustering_lcnstrs_merging(W, vertex_weights, CL, ML, normalized, verbosity)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    Y = 0;
    cnstr1 = [];
    cnstr2 = [];
    
    if ~isempty(CL)
        cnstr1 = CL(:,1);
        cnstr2 = CL(:,2);
    end
    
    temp = cell(1,1);
    temp{1,1} = ones(size(W,1),1);
    
    [W, vertex_weights, map, prev_clusters] = process_mls( W, vertex_weights, ML, temp );
    cnstr1 = map(cnstr1);
    cnstr2 = map(cnstr2);   
    
    [dML] = derive_mls_frm_cls( W, [cnstr1 cnstr2] );
    [W, vertex_weights, map_dml] = merge( W, vertex_weights, dML, temp );
    W = triu(W);
    W = W + W';
    cnstr1 = map_dml(cnstr1);
    cnstr2 = map_dml(cnstr2);
    
    [clusters,cuts,cheegers, vmin, W, vertex_weights, cnstr1, cnstr2] = ...
        spec_clustering_lcnstrs(W, vertex_weights, [cnstr1 cnstr2], [], normalized, verbosity);

    clusters = clusters(map_dml);
    clusters = clusters(map);

end