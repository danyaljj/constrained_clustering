function [comp,connected,sizes]=connectedComponents(W)
% Returns all connected components of the graph represented by weight
% matrix W.
%
% Usage: [comp,connected,sizes]=connectedComponents(W)
%
% (C)2010 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de

    l=1;
    [isCon, first_comp] = isConnected(W);
    sizes(l)=sum(first_comp);
    connected=isCon;

    comp=first_comp;
    while ~isCon
        ind=find(comp==0);
        Wpart=W(ind,ind);
        [isCon, first_comp] = isConnected(Wpart);
        l=l+1;
        sizes(l)=sum(first_comp);
        comp(ind)=l*first_comp;
    end
    
end
