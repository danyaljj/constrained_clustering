function Q = construct_cnstr_graph( cnstrs, n )
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    


    Q = sparse(n,n);
    if ~isempty(cnstrs)
        Q = sparse(cnstrs(:,1), cnstrs(:,2), 1, n,n) + sparse(cnstrs(:,2), cnstrs(:,1), 1, n,n);
    
        Q(Q>0) = 1;
    end

end