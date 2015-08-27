load tr11
kc = max(classid);
dat = unitnorm(logidf(dtm),2)';

for i = 1 : 10
	disp(sprintf('i = %d ...',i));

	tic; [pw,ci,ll,la] = kberns(dtm,kc);
	t(i,1) = toc;
	m(i,1) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = skberns(dtm,kc);
	t(i,2) = toc;
	m(i,2) = mi(cm(classid,ci));

	tic; [mob,ci,ll,la] = mixberns(dtm,kc);
	t(i,3) = toc;
	m(i,3) = mi(cm(classid,ci));

	tic; [mob,ci,ll,la,en] = daberns(dtm,kc);
	t(i,4) = toc;
	m(i,4) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = kmns(dtm,kc);
	t(i,5) = toc;
	m(i,5) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = skmns(dtm,kc);
	t(i,6) = toc;
	m(i,6) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = mixmns(dtm,kc);
	t(i,7) = toc;
	m(i,7) = mi(cm(classid,ci));

	tic; [mom,ci,ll,la,en] = damns(dtm,kc);
	t(i,8) = toc;
	m(i,8) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = kvmfs(dat,kc);
	t(i,9) = toc;
	m(i,9) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = skvmfs(dat,kc);
	t(i,10) = toc;
	m(i,10) = mi(cm(classid,ci));

	tic; [pw,ci,ll,la] = mixvmfs(dat,kc);
	t(i,11) = toc;
	m(i,11) = mi(cm(classid,ci));

	tic; [mov,ci,ll,la,en] = davmfs(dat,kc);
	t(i,12) = toc;
	m(i,12) = mi(cm(classid,ci));

end

