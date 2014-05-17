function hdp = hdp_setdata(hdp,dpindex,datass);

if nargin==2
  datass  = dpindex;
  dpindex = 1:hdp.numdp;
end

check_dpindex(hdp.numdp,dpindex);
if length(dpindex) ~= length(datass)
  error('dpindex and datass lengths do not match');
end

HELDOUT = 0;

for jj = 1:length(dpindex)
  if hdp.dpstate(dpindex(jj)) ~= HELDOUT
    error('Cannot set data for DPs that are not held out');
  end
  hdp.dp{dpindex(jj)}.numdata = size(datass{jj},2);
  hdp.dp{dpindex(jj)}.datass  = datass{jj};
end
