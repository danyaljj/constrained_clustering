function [idx centroids iter LCVQE time] = lcvqe(X, k, constraints,maxIter,centroids)
%LCVQE Modified Constrained vector quantization error (Pelleg & Baras, 2007)
% Modification on the k-means algorithm to penalize solutions that 
% violates pairwise constraints
% Implemented as described in the technical report except for the treatment of cannot link
% conList is a Cx3 matrix with the first and second values refering to 
% objects and the third column equal to 1 or -1, if the constraint
% is must-link or cannot-link
time = NaN;
tStart = tic;

numObjects = size(X,1);

if nargin == 3
	seeds = randsample( numObjects, k );
	centroids = X(seeds, :);
	maxIter=500;
end

distances = zeros( numObjects, k );
iter=1;
oldIdx = [];
idx = 1;

MLs = find( constraints(:,3) == 1 )';
CLs = find( constraints(:,3) == -1)';

while iter<=maxIter && ~isequal(oldIdx, idx) 
	%set of itens violating must and cannot link constraints
	GMLV = repmat(struct( 'item', [] ), k, 1);
	GCLV = repmat(struct( 'item', [] ), k, 1);

	oldIdx = idx;

	for c=1:k
		diffToCentroid = X - repmat(centroids(c,:), numObjects, 1);
		distances(:,c) = sum( diffToCentroid .^ 2, 2 ); 
	end

	[ minDists idx ] = min(distances, [], 2);

	
	%must link treatment
	for c=MLs
		s_1 = constraints(c, 1);
		s_2 = constraints(c, 2);
		c_j = idx(s_1);
		c_n = idx(s_2);

		%avoid useless calculations
		if c_j == c_n
			continue;
		end

		%check which choice minimize the error
		%keep in different cluster and violates the constraint
		case_a = 0.5*(distances(s_1, c_j) + distances(s_2, c_n)) + 0.25*(distances(s_1, c_n)+distances(s_2, c_j));
		%put s_2 in c_j
		case_b = 0.5*distances(s_1, c_j) + 0.5*distances(s_2, c_j);
		%put s_1 in c_n
		case_c = 0.5*distances(s_1, c_n) + 0.5*distances(s_2, c_n);

		cases = [ case_a case_b case_c ];
		[ mins idxMinCase ] = min(cases);
		
		switch(idxMinCase)

			case 1,%case a is minimal
				GMLV(c_j).item = [ GMLV(c_j).item; s_2 ];
				GMLV(c_n).item = [ GMLV(c_n).item; s_1 ];
				idx(s_1) = c_j;
				idx(s_2) = c_n;

			case 2,%case b is minimal
				idx(s_1) = c_j;
				idx(s_2) = c_j;

			case 3,%case c is minimal 
				idx(s_1) = c_n;
				idx(s_2) = c_n;

		end
	end
	
	if ~isempty(CLs)
		[ sortedDists idxSortedDists ] = sort( distances, 2 );
	end

	%cannot link treatment
	for c=CLs
		s_1 = constraints(c, 1);
		s_2 = constraints(c, 2);
		c_j = idx(s_1);
		c_n = idx(s_2);
		
		%avoid useless computations
		if c_j ~= c_n
			continue;
		end

		%setting the object farthest from the centroid of the cluster
		%and the neighbor cluster that is closest to this object
		r_j = s_1;
		MM_j = idxSortedDists(s_1, 2);
		closestObject = s_2;
		if( distances(s_2, c_j) > distances(s_1, c_j) )
			r_j = s_2;
			MM_j = idxSortedDists(s_2, 2);
			closestObject = s_1;
		end

		%
		%keep the violation
		case_a = 0.5*distances(s_1, c_j) + 0.5*distances(s_2, c_j) + 0.5*distances(r_j, MM_j);

		%test described by the author.. one of the cases is not checked (when assigning the object
		% closest to the centroid to another cluster, it seems counterintuitive but it can reduce the
		% lcvqe, my tests are below, and there it is covered
		%this case is for moving the farthest object from the c_j centroid to its second nearest cluster
		case_b = 0.5*distances( closestObject, c_j ) + 0.5*distances( r_j, MM_j);
		case_c = inf; %this is only used in my tests


		% my tests covering the case mentioned earlier, if used remember to CHANGE the switch below
		%move s_1 to the nearest neighbor cluster
		%case_b = 0.5*(distances(s_1, idxSortedDists(s_1,2)) + distances(s_2, c_n));
		%move s_2 to the nearest neighbor cluster
		%case_c = 0.5*(distances(s_1, c_j) + distances(s_2, idxSortedDists(s_2,2)));


		cases = [ case_a case_b case_c ];
		[ mins idxMinCase ] = min(cases);

		switch(idxMinCase)
			case 1,%case a is minimal				 
				GCLV(MM_j).item = [ GCLV(MM_j).item; r_j ];
				idx(s_1) = c_j;
				idx(s_2) = c_j;

			case 2,%case b is minimal
				%use this when using the authors test
				idx( closestObject ) = c_j;
				idx( r_j ) = MM_j;

				%use this when using my test
				%idx(s_1) = idxSortedDists(s_1,2);
				%idx(s_2) = c_n;

			case 3,%case c is minimal (only used in my tests)
				idx(s_1) = c_j;
				idx(s_2) = idxSortedDists(s_2,2);
		end

	end

	for c=1:k
		members=find(idx==c);
		coordsMembers = sum( X( members ,:) );
		coordsGMLV = sum( X(GMLV(c).item,:) );
		coordsGCLV = sum( X(GCLV(c).item,:) );
		n_j = length(members) + 0.5*length(GMLV(c).item) + length(GCLV(c).item);
		centroids(c,:) = (coordsMembers + 0.5*coordsGMLV + coordsGCLV) ./ n_j;
	end

	LCVQE = zeros(k,1);
	for c=1:k
		LCVQE(c) = 0.5*sum(distances(idx == c, c));
		sumML = 0;
		sumCL = 0;
		for itemViolated=GMLV(c).item'
			sumML = sumML + 0.5*distances(itemViolated, c);
		end
		for itemViolated=GCLV(c).item'
			sumCL = sumCL + distances(itemViolated, c);			
		end
		
		LCVQE(c) = LCVQE(c) + 0.5*sumML + 0.5*sumCL;
	end
	LCVQE = sum(LCVQE);

	time(iter) = toc(tStart);
	iter = iter + 1;
end


