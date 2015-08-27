
% test on tr11.mat
load tr11.mat
dat = unitnorm(logidf(dtm),2)';
kc = max(classid);
disp('test on tr11.mat ...');


[pw,ci,ll,la] = kberns(dtm,kc);
[pw,ci,ll,la] = kberns(dtm,kc,'ml');
[pw,ci,ll,la] = kberns(dtm,kc,'training','ml','maxi',10,'iprobs',pw);
[pw,ci,ll,la] = kberns(dtm,kc,'training','map','maxi',10,'ipart',classid);
pw0 = kberns(dtm,kc,'training','map','maxi',0,'ipart',classid);
[pw,ci] = kberns(dtm,kc,'training','map','maxi',1,'iprobs',pw0);
disp('kberns finished ...');

[pw,ci,ll,la] = skberns(dtm,kc);
[pw,ci,ll,la] = skberns(dtm,kc,'ml');
[pw,ci,ll,la] = skberns(dtm,kc,'training','map','maxi',10,'ipart',classid);
[pw,ci,ll,la] = skberns(dtm,kc,'training','ml','maxi',10,'iprobs',pw);
disp('skberns finished ...');

[mob,ci,ll,la] = mixberns(dtm,kc);
[mob,ci,ll,la] = mixberns(dtm,kc,'ml');
[mob,ci,ll,la] = mixberns(dtm,kc,'training','ml','maxi',10,'iparas',mob);
[mob,ci,ll,la] = mixberns(dtm,kc,'training','map','maxi',10,'ipart',classid);
mob0 = mixberns(dtm,kc,'map','maxi',0,'ipart',classid);
[mob,ci] = mixberns(dtm,kc,'map','maxi',1,'iparas',mob0);
disp('mixberns finished ...');

[pw,ci,ll,la] = bkberns(dtm,kc);
[pw,ci,ll,la] = bkberns(dtm,kc,'ml');
[pw,ci,ll,la] = bkberns(dtm,kc,'map','minn',20,'iprobs',pw);
[pw,ci,ll,la] = bkberns(dtm,kc,'training','map','minn',20,'ipart',classid);
disp('bkberns finished ...');


[pw,ci,ll,la] = kmns(dtm,kc);
[pw,ci,ll,la] = kmns(dtm,kc,'ml');
[pw,ci,ll,la] = kmns(dtm,kc,'map','maxi',10,'iprobs',pw);
[pw,ci,ll,la] = kmns(dtm,kc,'training','map','maxi',10,'ipart',classid);
pw0 = kmns(dtm,kc,'training','map','maxi',0,'ipart',classid);
[pw,ci] = kmns(dtm,kc,'training','map','maxi',1,'iprobs',pw0);
disp('kmns finished ...');

[pw,ci,ll,la] = skmns(dtm,kc);
[pw,ci,ll,la] = skmns(dtm,kc,'ml');
[pw,ci,ll,la] = skmns(dtm,kc,'training','map','maxi',10,'ipart',classid);
disp('skmns finished ...');

[mom,ci,ll,la] = mixmns(dtm,kc);
[mom,ci,ll,la] = mixmns(dtm,kc,'ml');
[mom,ci,ll,la] = mixmns(dtm,kc,'training','map','maxi',10,'iparas',mom);
[mom,ci,ll,la] = mixmns(dtm,kc,'training','map','maxi',10,'ipart',classid);
mom0 = mixmns(dtm,kc,'training','map','maxi',0,'ipart',classid);
[mom,ci] = mixmns(dtm,kc,'training','map','maxi',1,'iparas',mom0);
disp('mixmns finished ...');

[pw,ci,ll,la] = bkmns(dtm,kc);
[pw,ci,ll,la] = bkmns(dtm,kc,'ml');
[pw,ci,ll,la] = bkmns(dtm,kc,'training','map','minn',20,'iprobs',pw);
[pw,ci,ll,la] = bkmns(dtm,kc,'training','map','minn',20,'ipart',classid);
disp('bkmns finished ...');


[mu,ci,cs,ca] = kvmfs(dat,kc);
[mu,ci,cs,ca] = kvmfs(dat,kc,'maxi',10,'ipart',classid);
[mu,ci,cs,ca] = kvmfs(dat,kc,'maxi',10,'imeans',mu);
mu0 = kvmfs(dat,kc,'maxi',0,'ipart',classid);
[mu,ci] = kvmfs(dat,kc,'maxi',1,'imeans',mu0);
disp('kvmfs finished ...');

[mu,ci,cs,ca] = skvmfs(dat,kc);
[mu,ci,cs,ca] = skvmfs(dat,kc,'maxi',10,'ipart',classid);
[mu,ci,cs,ca] = skvmfs(dat,kc,'maxi',10,'imeans',mu);
disp('skvmfs finished ...');

[mov,ci,cs,ca] = mixvmfs(dat,kc);
[mov,ci,cs,ca] = mixvmfs(dat,kc,'maxi',10,'ipart',classid);
[mov,ci,cs,ca] = mixvmfs(dat,kc,'maxi',10,'iparas',mov);
mov0 = mixvmfs(dat,kc,'kappa',100,'maxi',0,'ipart',classid);
[mov,ci] = mixvmfs(dat,kc,'kappa',100,'maxi',1,'iparas',mov0);
disp('mixvmfs finished ...');

[mu,ci,cs,ca] = bkvmfs(dat,kc);
[mu,ci,cs,ca] = bkvmfs(dat,kc,'minn',20,'ipart',classid);
[mu,ci,cs,ca] = bkvmfs(dat,kc,'minn',20,'imeans',mu);
disp('bkvmfs finished ...');

