function hdp = hdp_setstate(hdp,hdpstate);

hdp.base.numclass = hdpstate.numclass;
hdp.base.classqq  = hdpstate.classqq;
hdp.dpstate       = hdpstate.dpstate;

for cp = 1:hdp.numconparam
  hdp.conparam{cp}.alpha = hdpstate.alpha(cp);
end
for jj = 1:hdp.numdp
  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);
  hdp.dp{jj}.datacc  = hdpstate.datacc{jj};
  hdp.dp{jj}.classnd = hdpstate.classnd{jj};
  hdp.dp{jj}.classnt = hdpstate.classnt{jj};
  hdp.dp{jj}.beta    = hdpstate.beta{jj};
  hdp.dp{jj}.alpha   = hdp.conparam{cp}.alpha;
  hdp.conparam{cp}.totalnd(tt) = sum(hdp.dp{jj}.classnd);
  hdp.conparam{cp}.totalnt(tt) = sum(hdp.dp{jj}.classnt);
end
