function hdp = dp_activate(hdp,dpindex,initcc);

check_dpindex(hdp.numdp,dpindex);

ACTIVE  = 2;
FROZEN  = 1;
HELDOUT = 0;

% initialize numclass and classqq
if ischar(initcc) & strcmp(initcc,'1perdp')
  oldclass = hdp.base.numclass;
  hdp = qq_addclass(hdp,length(dpindex));
elseif iscell(initcc)
  maxcc = max(cat(2,initcc{:}));
  if maxcc > hdp.base.numclass
    hdp = qq_addclass(hdp,maxcc - hdp.base.numclass);
  end
else
  if initcc > hdp.base.numclass
    hdp = qq_addclass(hdp,initcc - hdp.base.numclass);
  end
end

dpindex = sort(dpindex);
% initialize state of HDP
for kk = 1:length(dpindex)
  jj = dpindex(kk);
  pp = hdp.ppindex(jj);
  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);

  if hdp.dpstate(jj) == ACTIVE
    error('Trying to activate a DP that is already activated');
  end
  if pp > 0
    if hdp.dpstate(pp) ~= ACTIVE
      error('Ancestors of to be activated DPs has to be already activated');
    end
  end

  if hdp.dpstate(jj) == HELDOUT
    if ischar(initcc) & strcmp(initcc,'1perdp')
      hdp.dp{jj}.datacc = (oldclass+kk)*ones(1,hdp.dp{jj}.numdata);
    elseif iscell(initcc)
      hdp.dp{jj}.datacc = initcc{kk};
    else
      hdp.dp{jj}.datacc = ceil(rand(1,hdp.dp{jj}.numdata)*hdp.base.numclass);
    end
    hdp = qq_additems(hdp,hdp.dp{jj}.datacc,hdp.dp{jj}.datass);
    hdp.dp{jj}.classnd = full(bin(hdp.dp{jj}.datacc,1:hdp.base.numclass+1));
    hdp.dp{jj}.classnt = 0;
  end 
  if pp == 0
    hdp.dp{jj}.beta = ones(1,hdp.base.numclass+1)/(hdp.base.numclass+1);
  else
    hdp.dp{jj}.beta = hdp.dp{pp}.beta;
  end
  hdp.dp{jj}.alpha = hdp.conparam{cp}.alpha;
  hdp.dpstate(jj)  = ACTIVE;
end


for kk = length(dpindex):-1:1
  jj = dpindex(kk);
  pp = hdp.ppindex(jj);
  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);
  alpha = hdp.dp{jj}.alpha;

  if pp == 0
    hdp.dp{jj}.classnt = double(hdp.dp{jj}.classnd > 0);
  else
    hdp.dp{pp}.classnd = hdp.dp{pp}.classnd - hdp.dp{jj}.classnt;
    hdp.dp{jj}.classnt = randnumtable(alpha*hdp.dp{pp}.beta,hdp.dp{jj}.classnd);
    hdp.dp{pp}.classnd = hdp.dp{pp}.classnd + hdp.dp{jj}.classnt;
  end
  hdp.conparam{cp}.totalnd(tt) = sum(hdp.dp{jj}.classnd);
  hdp.conparam{cp}.totalnt(tt) = sum(hdp.dp{jj}.classnt);
end
 
