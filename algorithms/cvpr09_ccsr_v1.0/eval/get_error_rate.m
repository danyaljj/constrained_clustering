function [err corr] = get_error_rate(result, label)
%
% [ERR CORR]= GET_ERROR_RATE(RESULT, LABEL)
%
% To get the error rate for a clustering result, compared to the given
% clustering labels.
%
% ERR is the error rate of clustering and CORR is the nx2 vector which
% denote mapping the from the label of clustering result to given labels.


% check data integrity
l = length(result);
if l ~= length(label),
	error('- failed in data integrity check.');
end;	

result = reshape(result,[l 1]);
label = reshape(label,[l,1]);

result_unique = unique(result);
label_unique = unique(label);

% check the integrity of result
if length(result_unique) ~= length(label_unique),
        error('- The clustering result is not consistent with label.');
end;

n = length(result_unique);

% build the bipartite graph 
W = zeros(n);

for I = 1:n
	for J = 1:n
		W(I,J) = length(find(result==result_unique(I)&label==label_unique(J)));
		%W(J,I) = W(I,J);
	end;
end;

% find the maximum matching of the derived bipartite graph.
M = maximum_matching_bipartite(W);

idx = find(M>0);
[X Y] = ind2sub([n n],idx);

corr = [result_unique(X) label_unique(Y)];
err = 1-sum(W(idx))/l;

