function[centroid, pointsInCluster, assignment, clustersSize]= constrained_dpmeans(data, lambda, E, xi)
% function bpmeans()

% data
% lambda = 20;
% load tmp


% size(data,1) : number of the instances
% size(data,2) : dimensionality of the data

% data = X
% nbCluster = 5;
% load tmp

nbCluster = 1; % size(data,1)/10;
% nbClusterNow = 1;  % current number of clusters

% usage
% function[centroid, pointsInCluster, assignment]=
% myKmeans(data, nbCluster)
%
% Output:
% centroid: matrix in each row are the Coordinates of a centroid
% pointsInCluster: row vector with the nbDatapoints belonging to
% the centroid
% assignment: row Vector with clusterAssignment of the dataRows
%
% Input:
% data in rows
% nbCluster : nb of centroids to determine
%
% (c) by Christian Herta ( www.christianherta.de )
%
%
data_dim = length(data(1,:));
nbData  = length(data(:,1));

% init the centroids randomly
data_min = min(data);
data_max = max(data);
data_diff = data_max - data_min ;
% every row is a centroid
centroid = ones(nbCluster, data_dim) .* rand(nbCluster, data_dim);
for i=1 : 1 : length(centroid(:,1))
    centroid( i , : ) =   centroid( i , : )  .* data_diff;
    centroid( i , : ) =   centroid( i , : )  + data_min;
end
% end init centroids

% no stopping at start
pos_diff = 1.;

% main loop until
iterAll = 1;
while pos_diff > 0.0
    iterAll = iterAll + 1
    if iterAll > 20000
        disp('terminated by reaching the maximum number of iterations ')
        break;
    end
    %     nbClusterNow
    nbCluster
    % E-Step
    assignmentT = []; % per data
    % assign each datapoint to the closest centroid
    for d = 1 : length( data(:, 1) )
        min_diff = inf; 
        curAssignment = 0;
        
        for c = 1 : nbCluster;
            %             c
            diff2c = ( data( d, :) - centroid( c,:) );
            diff2c = sqrt( diff2c * diff2c' );
            if( min_diff >= diff2c)
                curAssignment = c;
                min_diff = diff2c;
            end
        end
        
%         if( min_diff < lambda )
%             % keep it ;
%         else
%             nbCluster =  nbCluster + 1;
%             disp([ 'Adding new clusters : min_diff '  num2str(min_diff) ] )
%             curAssignment = nbCluster;
%             centroid(end+1, :) = data( d, :); % (rand( 1, data_dim) .* data_diff) + data_min;
%             
%         end
        
        % assign the d-th dataPoint
        assignmentT = [ assignmentT; curAssignment];
    end

    Friends = zeros(length(data(:,1)), nbCluster); 
    Strangers = zeros(length(data(:,1)), nbCluster); 
    for d = 1:length(data(:,1))
        for d2 = 1:length(data(:,1))
            if E(d, d2) == 1
                Friends(d, assignmentT(d)) = Friends(d, assignmentT(d)) + 1;  
            elseif E(d, d2) == -1
                Strangers(d, assignmentT(d)) = Strangers(d, assignmentT(d)) + 1;  
            end 
        end
    end
    
    Friends
    Strangers 
    
    assignment = []; % per data
    % assign each datapoint to the closest centroid
    for d = 1 : length( data(:, 1) )
        min_diff = inf; 
        curAssignment = 0;
        
        for c = 1 : nbCluster;
            %             c
            diff2c = ( data( d, :) - centroid( c,:) );
            diff2c = sqrt( diff2c * diff2c' );
            
            if(  c <= size(Friends,2) ) 
                diff2c  = diff2c  - xi * ( Friends(d,c) - Strangers(d,c) ); 
            end 
            
            if( min_diff >= diff2c)
                curAssignment = c;
                min_diff = diff2c;
            end
        end
        
        if( min_diff < lambda )
            % keep it ;
        else
            nbCluster =  nbCluster + 1;
            disp([ 'Adding new clusters : min_diff '  num2str(min_diff) ] )
            curAssignment = nbCluster;
            centroid(end+1, :) = data( d, :); % (rand( 1, data_dim) .* data_diff) + data_min;
            
%             figure; hold on; 
%             for c = 1 : nbCluster-1;
%                 plot(centroid( c, 1), centroid( c, 2), 'xb'); 
%                 disp('plot another')
%             end
%             plot(centroid( c, 1), centroid( c, 2), 'xr'); 
%             pause; 
        end
        
        % assign the d-th dataPoint
        assignment = [ assignment; curAssignment];
    end

    
%     plotTheClusters(data, nbCluster,  centroid); 
%      pause; 
    
    
    % for the stoppingCriterion
    oldPositions = centroid;
    
    % M-Step
    % recalculate the positions of the centroids
    centroid = zeros(nbCluster, data_dim);
    pointsInCluster = zeros(nbCluster, 1);
    
    for d = 1: length(assignment)
        centroid( assignment(d),:) = data(d,:) + centroid( assignment(d),:)  ;
        pointsInCluster( assignment(d), 1 ) = pointsInCluster( assignment(d), 1 ) + 1;
    end
    
    for c = 1: nbCluster
        if( pointsInCluster(c, 1) ~= 0)
            try 
                centroid( c , : ) = centroid( c, : ) / pointsInCluster(c, 1);
            catch 
                disp('****************************************ERROR ***************************************')
                size(centroid)
                nbCluster
            end 
        else
            % set cluster randomly to new position
            %             centroid( c , : ) = (rand( 1, data_dim) .* data_diff) + data_min;
            centroid( c , : ) = [];
            nbCluster = nbCluster - 1; 
            disp('I need to remove this bitch! ')
            
        end
    end
    
%      plotTheClusters(data, nbCluster,  pointsInCluster); 
%      pause; 
    
    %stoppingCriterion
    if size(centroid,1 ) ~= size(oldPositions,1) || size(centroid,2) ~= size(oldPositions,2)
        pos_diff = 1; 
    else 
        pos_diff = sum (sum( (centroid - oldPositions).^2 ) );
    end 
    %     pos_diff
    clustersSize = nbCluster;
    if(pos_diff <= 0 )
        disp('terminated by reaching at the while threahold ')
    end 
    
end
