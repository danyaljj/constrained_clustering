function [ss, lmss]  = stirling(nn);

% returns the unsigned stirling numbers of the first kind.
% ss(nn,tt) for tt=1:nn

persistent maxnn allss logmaxss

if isempty(maxnn)
  maxnn = 1;
  allss = {1};
  logmaxss = 0;
end

if nn > maxnn
  allss{nn} = [];
  logmaxss(nn) = 0;
  for mm=maxnn+1:nn
    allss{mm} = [allss{mm-1}*(mm-1) 0] + [0 allss{mm-1}];
    mss = max(allss{mm});
    allss{mm} = allss{mm}/mss;
    logmaxss(mm) = logmaxss(mm-1) + log(mss);
  end
  maxnn = nn;
end

ss = allss{nn};
lmss = logmaxss(nn);
