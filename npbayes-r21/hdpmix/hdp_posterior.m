function [sample,hdp,lik]=hdp_posterior(hdp,numburnin,numsample,numspace,...
                                        doconparam,dolik,dodebug,fid);
% runs hdp, gets samples for classqq, totalnt.  These are enough to form
% predictive likelihoods for unseen data.

totiter = numburnin + numsample*numspace;
sample  = cell(1,numsample);
lik     = zeros(1,numburnin+numsample*numspace);

countdown(fid,'',0,numburnin + numsample*numspace);
[hdp ll] = hdp_iterate(hdp,numburnin,doconparam,dolik,dodebug);
lik(1:numburnin) = ll;
countdown(fid,'\nburn in: ',numburnin);
for samp = 1:numsample
  [hdp ll] = hdp_iterate(hdp,numspace,doconparam,dolik,dodebug);
  lik(numburnin+(samp-1)*numspace+(1:numspace)) = ll;
  countdown(fid,['\npostsample ' num2str(samp) ': '],numspace);
  sample{samp} = hdp_getstate(hdp);
end

if dolik == 0
  clear lik
end
