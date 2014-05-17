function [ dist ] = gaussianDifference( x,y )

dist = x - y;
dist = sum(dist*dist');
dist = sqrt(dist);

end

