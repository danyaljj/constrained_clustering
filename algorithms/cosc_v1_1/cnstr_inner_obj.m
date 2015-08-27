function obj = cnstr_inner_obj( D, FctVal, r2, vec, gamma, volQ, g)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    obj = sum(abs(D*g))+gamma*volQ*(max(g)-min(g))/2 - g'*r2/2 - FctVal*g'*vec/2;
    
end