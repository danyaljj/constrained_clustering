function [cc, numclass] = randcrp(alpha, numdata);
% generates a CRP partition of numdata items, with concentration parameter
% alpha.

cc = zeros(1,numdata);
weights  = alpha;
numclass = 0;

for ii = 1:numdata
  cc(ii) = randmult(weights);
  if cc(ii) > numclass
    weights = [weights(1:numclass) 1 alpha];
    numclass = cc(ii);
  else
    weights(cc(ii)) = weights(cc(ii)) + 1;
  end
end
