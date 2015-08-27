function ccut = compute_cheeger_cut(W, gdeg, Y)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    labels = unique(Y);
    assert( length(labels) == 2, 'only binary partitions are allowed');
    
    ix = sum(Y == labels(1) );
    ixc = sum(Y == labels(2) );
    cut = sum(sum(W(ix, ixc)));
    vol = sum(gdeg(ix));
    volc = sum(gdeg(ixc));
    
    ccut = cut/ min(vol, volc);

end