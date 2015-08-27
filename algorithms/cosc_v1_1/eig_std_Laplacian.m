function [v2, l2, ones_C] = eig_std_Laplacian(W,normalized,deg)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    W = triu(W);
    W = W+W';
    
    D=spdiags(sum(W,2),0,size(W,1),size(W,1));
    opts.disp=0;
    options.tol = 1E-6;
    options.maxit=20;
    options.issym = 1;
    if (normalized)
        [eigvec,eigval]= eigs(D-W, spdiags(deg,0,size(W,1),size(W,1)),2,'SA',opts);
        %[eigvec,eigval]= eigs(D-W, diag(deg),2,'SA',opts);
    else
        [eigvec,eigval]= eigs(D-W, 2,'SA',opts);
    end
    v2 = eigvec(:,2);
    l2 = eigval(2,2);
    
    start = createClustersGeneral(v2,W,normalized,-1,2,deg,true);
    if (sum(start)>sum(start==0)) start=1-start; end
    start=start/sum(start);
    ones_C = start;
    
end
