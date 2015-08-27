function[centroid, pointsInCluster, assignment, clustersSize]= constrained_kmeans(data, E, k)

[num dim] = size(data);

pointsInCluster = []; 

% init the centroids randomly
data_min = min(data);
data_max = max(data);
data_diff = data_max - data_min;

% every row is a centroid
centroid = rand(k, dim);
for i=1 : 1 : length(centroid(:,1))
    centroid(i,:) = centroid(i,:) .* data_diff;
    centroid(i,:) = centroid(i,:) + data_min;
end

% no stopping at start
pos_diff = 1.;

% main loop until
iterAll = 0;
while pos_diff > 0.0
    iterAll = iterAll + 1
    if iterAll >= 20000
        disp('terminated by reaching the maximum number of iterations ')
        break;
    end
    
    % E-Step
    assignment = []; % per data
    % assign each datapoint to the closest centroid
    for d = 1:1:num
        min_diff = inf;
        curAssignment = 0;
        for c = 1:1:k
            diff2c = (data(d,:) - centroid(c,:));
            diff2c = sqrt(diff2c*diff2c');
            if( min_diff >= diff2c)
                %TODO check constraints
                violation = 0;
                for dd=1:1:length(assignment)
                    if(E(d,dd)<0 || E(dd,d)<0)
                        %check if a point has been assigned to the cluster
                        %that violated a NOT constraint
                        if(assignment(dd)==c)
                            violation = 1;
                        end
                    end
                    if(E(d,dd)>0 || E(dd,d) > 0)
                        %check if a MUST constraint has been assigned to a different cluster.
                        if(~(assignment(dd)==c))
                            violation = 1;
                        end
                    end
                end
                
                if(violation==0)
                    %make assignment
                    curAssignment = c;
                    min_diff = diff2c;
                end
            end
        end
        
        if curAssignment == 0
           disp('impossible to assign the points!!!') 
           return
        end
        assignment = [ assignment; curAssignment];
    end
    
    assignment
    
    % for the stoppingCriterion
    oldPositions = centroid;
    
    % M-Step
    % recalculate the positions of the centroids
    centroid = zeros(dim, k);
    pointsInCluster = zeros(k, 1);
    
    for d = 1:1:length(assignment)
        centroid( assignment(d),:) = data(d,:) + centroid( assignment(d),:);
        pointsInCluster(assignment(d),1) = pointsInCluster(assignment(d),1) + 1;
    end
    
    for c = 1:1:k;
        if( pointsInCluster(c, 1) ~= 0)
            centroid( c , : ) = centroid( c, : ) / pointsInCluster(c, 1);
        else
            % set cluster randomly to new position
            centroid( c , : ) = (rand( 1, dim) .* data_diff) + data_min;
        end
    end
    
    %stoppingCriterion
    if size(centroid,1) ~= size(oldPositions,1) || size(centroid,2) ~= size(oldPositions,2)
        pos_diff = 1;
    else
        pos_diff = sum(sum((centroid-oldPositions).^2));
    end
    
    if(pos_diff <= 0 )
        disp('terminated by reaching at the while threahold ')
    end
end