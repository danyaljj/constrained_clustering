function [ dist ] = multDifference( x,y )
 
dist = 0;
%penalty = 10;

for i = 1:1:length(x)
   num = x(i);
   den = y(i);
   frac = num / den;
   %if y(i) == 0
   %    if x(i) ~= 0
           %dist = inf;
           %break;
   %        dist = dist + penalty;
   %        continue;
   %    else
   %        continue;
   %    end
   %end
   
   %if x(i) == 0
   %   if y(i) ~= 0
           %dist = inf;
           %break;
   %        dist = dist + penalty;
   %        continue;
   %   end
   %end

   dist = dist + x(i)*log(frac);

end

end

