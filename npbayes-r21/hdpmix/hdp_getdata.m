function datass = hdp_getdata(hdp,dpindex);

if nargin==1
  dpindex = 1:hdp.numdp;
end

check_dpindex(hdp.numdp,dpindex);

datass  = cell(1,length(dpindex));

for jj = 1:length(dpindex)
  datass{jj}  = hdp.dp{dpindex(jj)}.datass;
end

