function [hdp,dpindex] = hdp_adddp(hdp,numdp,pp,cp);

HELDOUT = 0;

dpindex   = hdp.numdp + (1:numdp);
hdp.numdp = hdp.numdp + numdp;
for ii = 1:numdp
  jj = dpindex(ii);
  tt = hdp.conparam{cp}.numdp + 1;
  hdp.dpstate(jj)              = HELDOUT;
  hdp.ppindex(jj)              = pp;
  hdp.cpindex(jj)              = cp;
  hdp.ttindex(jj)              = tt;
  hdp.conparam{cp}.numdp       = tt;
  hdp.conparam{cp}.totalnd(tt) = 0;
  hdp.conparam{cp}.totalnt(tt) = 0;
  hdp.dp{jj}.datacc            = [];
  hdp.dp{jj}.classnd           = 0;
  hdp.dp{jj}.classnt           = 0;
  hdp.dp{jj}.beta              = 1;
  hdp.dp{jj}.alpha             = [];
  hdp.dp{jj}.numdata           = 0;
  hdp.dp{jj}.datass            = [];

end
