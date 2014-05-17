function lik = hdp_predict(hdp,postsample,trainindex,testindex,...
                           numburnin,numsample,doconparam,dodebug,fid);

numgroup      = length(testindex);
numpostsample = length(postsample);
lik           = zeros(numgroup,numpostsample,numsample);
countdown(fid,'',0,numgroup*numpostsample*(numburnin+numsample));

for ps = 1:numpostsample
  hdp = hdp_setstate(hdp,postsample{ps});
  hdp = dp_freeze(hdp,trainindex);
  hdp = dp_activate(hdp,testindex,'1perdp');
  [hdp ll] = feval(hdp.func.predict,hdp,numburnin,numsample,doconparam,...
      testindex,dodebug);
  preamble = sprintf('samples %d: ',ps);
  countdown(fid,['\n' preamble],numgroup*(numsample+numburnin));
  lik(:,ps,:) = reshape(ll',[numgroup 1 numsample]);
end
