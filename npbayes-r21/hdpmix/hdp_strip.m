function hdp = hdp_strip(hdp)

HELDOUT = 0;

hdp.dpstate = HELDOUT*ones(1,hdp.numdp);
for jj = 1:hdp.numdp
  hdp.dp{jj}.datacc  = [];
  hdp.dp{jj}.classnd = [];
  hdp.dp{jj}.classnt = [];
  hdp.dp{jj}.beta    = [];
  hdp.dp{jj}.numdata = 0;
  hdp.dp{jj}.datass  = [];
end
for cp = 1:hdp.numconparam
  hdp.conparam{cp}.alpha   = hdp.conparam{cp}.alphaa/hdp.conparam{cp}.alphab;
  hdp.conparam{cp}.totalnd = zeros(1,hdp.conparam{cp}.numdp);
  hdp.conparam{cp}.totalnt = zeros(1,hdp.conparam{cp}.numdp);
end

