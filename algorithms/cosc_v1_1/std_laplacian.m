function L = std_laplacian(W)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
    D=spdiags(sum(W,2),0,size(W,1),size(W,1));
    L = D-W;

end