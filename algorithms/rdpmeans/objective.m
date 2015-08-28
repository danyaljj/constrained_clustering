function [ obj ] = objective(data, assignment, centroid, E, xi, lambda, distancef)

obj = 0;

if(nargin < 6)
    distancef ='Gaussian';
end

for d=1 : length(data(:,1))
   
   %obj = obj + sqrt(sum((data(d,:)-centroid(assignment(d),:)).^2));
   
   obj = 0;
   if(strcmp(distancef,'Gaussian'))
        obj = gaussianDifference(data(d,:), centroid(assignment(d),:));
   else
        obj = multDifference(data(d,:), centroid(assignment(d),:));
   end
   
   ptsInCluster = find(assignment == assignment(d));
   
   friends = 0;
   strangers = 0;
   for p =1: length(ptsInCluster)
       if E(ptsInCluster(p),d) == 1
           friends = friends + 1;
       elseif E(ptsInCluster(p),d) == -1
           strangers = strangers + 1;
       end
   end
   
   obj = obj - xi*(friends - strangers);
   %disp(['Friends = '  num2str(friends) ]); 
   %disp(['Strangers = '  num2str(strangers) ]);
   %disp(['obj = '  num2str(obj) ]);
   
end

obj = obj + lambda*size(centroid,1);

end

