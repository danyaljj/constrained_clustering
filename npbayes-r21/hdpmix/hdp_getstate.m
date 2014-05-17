function hdpstate = hdp_getstate(hdp);

hdpstate.numclass = hdp.base.numclass;
hdpstate.classqq  = hdp.base.classqq;
hdpstate.dpstate  = hdp.dpstate;
hdpstate.datacc   = cell(1,hdp.numdp);
hdpstate.classnd  = cell(1,hdp.numdp);
hdpstate.classnt  = cell(1,hdp.numdp);
hdpstate.beta     = cell(1,hdp.numdp);
hdpstate.alpha    = zeros(1,hdp.numconparam);

% collect posterior statistics stuff
for jj = 1:hdp.numdp
  hdpstate.datacc{jj}  = hdp.dp{jj}.datacc;
  hdpstate.classnd{jj} = hdp.dp{jj}.classnd;
  hdpstate.classnt{jj} = hdp.dp{jj}.classnt;
  hdpstate.beta{jj}    = hdp.dp{jj}.beta;
end 
for cp = 1:hdp.numconparam
  hdpstate.alpha(cp)   = hdp.conparam{cp}.alpha;
end


