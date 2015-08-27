function bcut = bal_cut( W, vertex_weights, Y)
% Usage: bcut = bal_cut( W, vertex_weights, Y)
% Computes balanced cut where the balance is specified by vertex_weights
% and partition by Y. Also works for multiple partitions.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    classes = unique(Y);
    k = length(classes);
    
    bcut = 0;
    n = size(W,1);
    
    for i=1:k
    
        jx = find(Y==classes(i));
        jxc = setdiff(1:n, jx);
        
        cut = sum(sum( W(jx, jxc) ));
        vol = sum(vertex_weights(jx));        

        bcut = bcut + cut/vol;
        
    end

end