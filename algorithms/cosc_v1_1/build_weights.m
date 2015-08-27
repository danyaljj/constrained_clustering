% Builds a KNN-Graph or epsilon-Neighbourhood Graph
function K=build_weights(points,graphType,adaptive,numKNN, scale)
% Usage: K=buildWeights(points,graphType,adaptive,numKNN)
%
% Input: 
% points        - the coordinates of points in Euclidean space. It is a num x dim matrix.
% graphType     - <2: symmetric KNN; otherwise epilson neighborhood
% adaptive      - true 
% numKNN        - number of neighbors to connect to
% scale         - 4 in the paper
% 
% (C)2010 Thomas Buehler, Syama Sundar Rangapuram, Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de


    % epsilon for graphType >= 2
    eps=0.3;
    
    
    %numKNN=10;
    if nargin<4
		numKNN=10;
    end

    % Compute the squares of the pairwise distances. dist2 is num x num
    % matrix
    dist2 =DistEuclideanPiotrDollar(points,points); % squared distances
    
    num =size(points,1);
    
    if(graphType<2)
        
        % ith column of SD contains neighbors sorted according to their
        % distances to ith point
        [SD,IX]=sort(dist2,1);
        KNN     = IX(2:numKNN+1,:);
        
        % KNNDist is numKNN x n matrix containing the (squares of) distances to numKNN neighbors (excluding self) 
        KNNDist = SD(2:numKNN+1,:);
        
        if (~adaptive)
            gamma_squared=mean(mean(KNNDist))*ones(1,num);
        else
            % find the square of the distance of the current point to its
            % numKNNth neighbor
            %gamma_squared=mean(KNNDist);%*ones(1,num);%KNNDist(numKNN,:);%mean(KNNDist);
            gamma_squared=KNNDist(numKNN,:);%mean(KNNDist);
        end
            
        % get kNN weight matrix
        K = sparse(num,num);
        for i=1:num
            K(KNN(:,i),i)=exp(-scale/(gamma_squared(i))*KNNDist(:,i));
        end
        % note that K is not symmetric yet , now we symmetrize K 
        if(graphType==1) K=(K+K')+abs(K-K'); K=0.5*K; end
        if(graphType==0) K=(K+K')-abs(K-K'); K=0.5*K; end 
        
        K=K-spdiags(diag(K),0,num,num);
    else
        K = exp(-1/(2*gamma^2)*dist2).*(dist2 < eps^2 & dist2~=0);
    end

end