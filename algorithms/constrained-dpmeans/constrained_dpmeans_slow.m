function[centroid, pointsInCluster, assignment, clustersSize, ... 
    objs, pointsAll, centroindsAll]= constrained_dpmeans_slow(data, lambda, E, rate, xi0, distancef)

if(nargin < 6)
    distancef ='Gaussian';
end

nbCluster = 1; 
data_dim = length(data(1,:));
nbData  = length(data(:,1));
objs = [];

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
pointsAll = {}; 
centroindsAll = {}; 
xi = 0; 

assignment = ones(size(data,1),1); % per data
assignment_old = assignment;
% main loop until
iterAll = 0;
while pos_diff > 0.0
    iterAll = iterAll + 1
    if iterAll > 20000
        disp('terminated by reaching the maximum number of iterations ')
        break;
    end

    [iterAll nbCluster]

    for d = randperm(length( data(:, 1) ))
        min_diff = inf; 
        curAssignment = 0;
        currd = inf;
        
        for c = 1 : nbCluster;
            % finding friends and strangers 
            friends = 0;
            strangers = 0; 
            
            if(iterAll > 1)
                ptsInCluster = find(assignment == c);
                %ptsInCluster = 1:length(assignment); 
				friends = 0;
                strangers = 0;
                for p =1: length(ptsInCluster)
                    if E(ptsInCluster(p),d) || E(d,ptsInCluster(p)) == 1
                        friends = friends + 1;
                    elseif E(ptsInCluster(p),d) || E(d,ptsInCluster(p)) == -1
                        strangers = strangers + 1;
                    end
                end
            end

            diff2c = 0;
            if(strcmp(distancef,'Gaussian'))
                diff2c = gaussianDifference(data(d,:), centroid(c,:));
            else
                diff2c = multDifference(data(d,:), centroid(c,:));
            end
            
            diff2c = diff2c  - xi * (friends - strangers); 
            
            if iterAll > 1
                if c==assignment(d)
                    currd = diff2c;
                end
            end
            
            if( min_diff >= diff2c)
                curAssignment = c;
                min_diff = diff2c;
            end
        end
        
        if min_diff > currd
            dis(['increase in assignment' num2str(min_diff) num2str(currd)]);
        end
        
        if( min_diff > lambda )
            nbCluster = nbCluster + 1;
            disp([ 'Adding new clusters : min_diff '  num2str(min_diff) ] );
            curAssignment = nbCluster;
            centroid(end+1, :) = data( d, :); % (rand( 1, data_dim) .* data_diff) + data_min
        end
        
        % assign the d-th dataPoint
        assignment(d) = curAssignment;
        
        if iterAll>1
        %calculate objective
            objs = [objs objective(data, assignment, centroid, E, xi, lambda, distancef)];
            pointsAll {end+1} = assignment; 
            centroindsAll{end+1} = centroid;
            
            %if length(objs) > 1 && objs(end-1) < objs(end)
            %    disp(['objective increased here!!!!!'  num2str(currd) num2str(min_diff) ])
            %end
        end
    end
    
    disp(['Assignment of points finished!! '  num2str(length(objs)) ])
    
    % for the stoppingCriterion
    oldPositions = centroid;
    assignment_old = assignment;
    
    % M-Step
    % recalculate the positions of the centroids
    centroid = zeros(nbCluster, data_dim);
    pointsInCluster = zeros(nbCluster, 1);
    
    for d = 1: length(assignment)
        centroid( assignment(d),:) = data(d,:) + centroid( assignment(d),:)  ;
        pointsInCluster( assignment(d), 1 ) = pointsInCluster( assignment(d), 1 ) + 1;
    end
    
    add = 0; 
    e = nbCluster; 
    for c = 1: e
        if( pointsInCluster(c, 1) ~= 0)
            try 
                centroid( c - add , : ) = centroid( c - add, : ) / pointsInCluster(c, 1);
            catch 
                disp('****************************************ERROR ***************************************')
            end 
        else
            % set cluster randomly to new position
            %             centroid( c , : ) = (rand( 1, data_dim) .* data_diff) + data_min;
            centroid( c - add , : ) = [];
            nbCluster = nbCluster - 1; 
            disp('Cluster removed!')
            add = add + 1;
            ind = find(assignment > c ); 
            assignment(ind) = assignment(ind) - 1; 
        end
    end
    
    %stoppingCriterion
    if size(centroid,1 ) ~= size(oldPositions,1) || size(centroid,2) ~= size(oldPositions,2)
        pos_diff = 1; 
    else 
        pos_diff = sum (sum( (centroid - oldPositions).^2 ) );
    end 

    clustersSize = nbCluster;
    if(pos_diff <= 0 )
        disp('terminated by reaching at the while threshold ')
    end 
    xi = xi0 * (rate)^(iterAll); 
    disp(['Current value of xi = ' num2str(xi) ])
end
