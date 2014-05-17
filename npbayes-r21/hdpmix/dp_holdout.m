function hdp = dp_holdout(hdp,dpindex);

check_dpindex(hdp.numdp,dpindex);

HELDOUT = 0;

for jj = -sort(-dpindex)
  if hdp.dpstate(jj) == HELDOUT
    error('DP already heldout');
  end

  pp = hdp.ppindex(jj);
  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);

  hdp = qq_delitems(hdp,hdp.dp{jj}.datacc,hdp.dp{jj}.datass);
  if pp > 0
    hdp.dp{pp}.classnd = hdp.dp{pp}.classnd - hdp.dp{jj}.classnt;
  end
  hdp.dp{jj}.datacc  = [];
  hdp.dp{jj}.classnd = [];
  hdp.dp{jj}.classnt = [];
  hdp.dp{jj}.beta    = [];
  hdp.dp{jj}.alpha   = [];
  hdp.conparam{cp}.totalnd(tt) = 0;
  hdp.conparam{cp}.totalnt(tt) = 0;
  hdp.dpstate(jj)    = HELDOUT;
end

check_dpstate(hdp.ppindex,hdp.dpstate);
