function [cut,cheeger,cutParts]= computeMultiCut(W,allClusters,normalized)
% Evaluates the multicut versions of Ratio/Normalized Cut and
% Ratio/Normalized Cheeger Cut
%
% Usage: [cut,cheeger,cutParts]= computeMultiCut(W,allClusters,normalized)
%
% (C)2010-11 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de

	if(normalized)
		deg=sum(W,2);
	else
		deg=ones(size(W,1),1);
	end
    
	labels=unique(allClusters);
    cut=0;
	cheeger=0;
    for k=1:length(labels)

        clustersM=zeros(size(allClusters,1),1);
        clustersM(allClusters==labels(k))=1;
        cutParts(k) = computeCutValue(clustersM,W,normalized,deg);
        
    end
    cut=sum(cutParts);
    cheeger=max(cutParts);
    
end