function hdp = hdp_randconparam(hdp,numiter);

HELDOUT = 2;
for cp = 1:hdp.numconparam
  hdp.conparam{cp}.alpha = randconparam(hdp.conparam{cp}.alpha,...
      hdp.conparam{cp}.totalnd,hdp.conparam{cp}.totalnt,...
      hdp.conparam{cp}.gammaa,hdp.conparam{cp}.gammab,numiter);
end

for jj = 1:hdp.numdp
  cp = hdp.cpindex(jj);
  if hdp.dpstate(jj)~=HELDOUT
    hdp.dp{jj}.alpha = hdp.conparam{cp}.alpha;
  end
end
