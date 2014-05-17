function hdp = hdp_deldp(hdp,dpindex);

check_dpindex(hdp.numdp,dpindex);

HELDOUT = 0;

for jj = 1:hdp.numdp
  pp = hdp.ppindex(jj);
  if ~isempty(find(pp==dpindex)) & isempty(find(jj==dpindex))
    error('Need to delete all descendents of DPs as well');
  end
end

for jj = -sort(-dpindex)
  if hdp.dpstate(jj) ~= HELDOUT
    error('DPs to be deleted must first be held out');
  end

  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);

  hdp.conparam{cp}.numdp = hdp.conparam{cp}.numdp - 1;
  hdp.conparam{cp}.totalnd(tt) = [];
  hdp.conparam{cp}.totalnt(tt) = [];
  for j2 = 1:hdp.numdp
    if hdp.cpindex(j2) == cp & hdp.ttindex(j2) > tt
      hdp.ttindex(j2) = hdp.ttindex(j2) - 1;
    end
  end
  hdp.dp(jj)      = [];
  hdp.ppindex(jj) = [];
  hdp.cpindex(jj) = []
  hdp.ttindex(jj) = []
  hdp.dpstate(jj) = []
  hdp.numdp       = hdp.numdp - 1;
end

