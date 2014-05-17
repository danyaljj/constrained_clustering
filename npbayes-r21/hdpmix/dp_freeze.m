function hdp = hdp_setfrozen(hdp,dpindex);

check_dpindex(hdp.numdp,dpindex);

ACTIVE = 2;
FROZEN = 1;

for jj = dpindex
  if hdp.dpstate(jj)~=ACTIVE
    error('Can only freeze a DP that is activated and not frozen');
  end
  hdp.dp{jj}.alpha = [];
  hdp.dp{jj}.beta  = [];
  hdp.dpstate(jj)  = FROZEN;
end

check_dpstate(hdp.ppindex,hdp.dpstate);
