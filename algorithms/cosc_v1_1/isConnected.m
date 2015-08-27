function [connected,components]=isConnected(W)
% Checks whether a graph is connected.
%
% Usage: [connected,components]=isConnected(W)
%
% (C)2010 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de

	A = W>0; % adjacency matrix

	alreadyseen = zeros(size(W,1),1);

	currentCandidates=1;

	while ~isempty(currentCandidates)
		candidates= (sum(A(:,currentCandidates),2)>0);
		alreadyseen(currentCandidates)=1;
		currentCandidates=find(candidates-alreadyseen>0);
	end

	connected = sum(alreadyseen)==size(W,2);
    
    components=alreadyseen;
    
end
