function hdp = hdp_init(func,ppindex,cpindex,hh,alphaa,alphab);

check_ppindex(length(ppindex),ppindex);
check_cpindex(max(cpindex),cpindex);

HELDOUT = 0;

hdp.func          = func;
hdp.numdp         = length(ppindex);
hdp.numconparam   = length(alphaa);
hdp.base.hh       = hh;
hdp.base.classqq  = feval(hdp.func.newclass,hdp.base.hh);
hdp.base.numclass = 0;
hdp.conparam      = cell(1,hdp.numconparam);
hdp.dp            = cell(1,hdp.numdp);
hdp.dpstate       = HELDOUT*ones(1,hdp.numdp);
hdp.ppindex       = ppindex;
hdp.cpindex       = cpindex;
hdp.ttindex       = zeros(1,hdp.numdp);

tt = zeros(1,hdp.numconparam);
for jj = 1:hdp.numdp
  cp              = cpindex(jj);
  tt(cp)          = tt(cp) + 1;
  hdp.ttindex(jj) = tt(cp);
end

for jj = 1:hdp.numdp
  hdp.dp{jj}.datacc  = [];
  hdp.dp{jj}.classnd = 0;
  hdp.dp{jj}.classnt = 0;
  hdp.dp{jj}.beta    = 1;
  hdp.dp{jj}.alpha   = [];
  hdp.dp{jj}.numdata = 0;
  hdp.dp{jj}.datass  = [];
end
for cp = 1:hdp.numconparam
  hdp.conparam{cp}.alphaa  = alphaa(cp);
  hdp.conparam{cp}.alphab  = alphab(cp);
  hdp.conparam{cp}.numdp   = sum(cpindex==cp);
  hdp.conparam{cp}.alpha   = hdp.conparam{cp}.alphaa/hdp.conparam{cp}.alphab;
  hdp.conparam{cp}.totalnd = zeros(1,hdp.conparam{cp}.numdp);
  hdp.conparam{cp}.totalnt = zeros(1,hdp.conparam{cp}.numdp);
end


