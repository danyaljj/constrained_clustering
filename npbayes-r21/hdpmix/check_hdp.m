function check_hdp(hdp);

ACTIVE = 2;
FROZEN = 1;
HELDOUT = 0; 

% Check for lengths
if length(hdp.conparam)~=hdp.numconparam
  error('length(hdp.conparam)~=hdp.numconparam');
end
if length(hdp.dp)~=hdp.numdp
  error('length(hdp.dp)~=hdp.numdp');
end
if length(hdp.dpstate)~=hdp.numdp
  error('length(hdp.dpstate)~=hdp.numdp');
end
if length(hdp.ppindex)~=hdp.numdp
  error('length(hdp.ppindex)~=hdp.numdp');
end
if length(hdp.cpindex)~=hdp.numdp
  error('length(hdp.cpindex)~=hdp.numdp');
end
if length(hdp.ttindex)~=hdp.numdp
  error('length(hdp.ttindex)~=hdp.numdp');
end

% Check indices
check_ppindex(hdp.numdp,hdp.ppindex);
check_cpindex(hdp.numconparam,hdp.cpindex);
check_ttindex(hdp.cpindex,hdp.ttindex);
if max(hdp.cpindex)~=hdp.numconparam
  error('max(hdp.cpindex)!=hdp.numconparam');
end

% Check conparams
for kk = 1:hdp.numconparam
  if length(hdp.conparam{kk}.alphaa) < 0
    error(sprintf('length(hdp.conparam{%d}.alphaa)<0',kk));
  end
  if length(hdp.conparam{kk}.alphab) < 0
    error(sprintf('length(hdp.conparam{%d}.alphab)<0',kk));
  end
  if length(hdp.conparam{kk}.alpha) < 0
    error(sprintf('length(hdp.conparam{%d}.alpha)<0',kk));
  end
  if length(hdp.conparam{kk}.totalnd)~=hdp.conparam{kk}.numdp
    error(sprintf('length(hdp.conparam{%d}.totalnd)!=hdp.conparam{%d}.numdp',...
                  kk,kk));
  end
  if length(hdp.conparam{kk}.totalnt)~=hdp.conparam{kk}.numdp
    error(sprintf('length(hdp.conparam{%d}.totalnt)!=hdp.conparam{%d}.numdp',...
                  kk,kk));
  end
end
 
% check DPs
check_dpstate(hdp.ppindex,hdp.dpstate);
for jj = hdp.numdp:-1:1
  pp = hdp.ppindex(jj);
  cp = hdp.cpindex(jj);
  tt = hdp.ttindex(jj);
  if size(hdp.dp{jj}.datass,2)~=hdp.dp{jj}.numdata
    error(sprintf('size(hdp.dp{%d}.datass,2)!=hdp.dp{%d}.numdata',jj,jj));
  end

  if hdp.dpstate(jj) == FROZEN | hdp.dpstate(jj) == ACTIVE
    if length(hdp.dp{jj}.classnd)~=hdp.base.numclass+1
      error(sprintf('length(hdp.dp{%d}.classnd)!=hdp.base.numclass+1',jj));
    end
    if length(hdp.dp{jj}.classnt)~=hdp.base.numclass+1
      error(sprintf('length(hdp.dp{%d}.classnt)!=hdp.base.numclass+1',jj));
    end
    if length(hdp.dp{jj}.datacc)~=hdp.dp{jj}.numdata
      error(sprintf('length(hdp.dp{%d}.datacc)!=hdp.dp{%d}.numdata',jj,jj));
    end
    if any(hdp.dp{jj}.classnt<0)
      error(sprintf('hdp.dp{%d}.classnt<0',jj));
    end
    if any(hdp.dp{jj}.classnd<hdp.dp{jj}.classnt)
      error(sprintf('hdp.dp{%d}.classnd<hdp.dp{%d}.classnt',jj,jj));
    end
    if hdp.conparam{cp}.totalnd(tt)~=sum(hdp.dp{jj}.classnd)
      error(sprintf('sum(hdp.dp{%d}.classnd)~=totalnd',jj));
    end
    if hdp.conparam{cp}.totalnt(tt)~=sum(hdp.dp{jj}.classnt)
      error(sprintf('sum(hdp.dp{%d}.classnt)~=totalnt',jj));
    end
    hdp = qq_delitems(hdp,hdp.dp{jj}.datacc,hdp.dp{jj}.datass);
  else
    if hdp.conparam{cp}.totalnd(tt)~=0
      error(sprintf('hdp.conparam{%d}.totalnd{%d}!=0',cp,tt));
    end
    if hdp.conparam{cp}.totalnt(tt)~=0
      error(sprintf('hdp.conparam{%d}.totalnt{%d}!=0',cp,tt));
    end
  end
  if hdp.dpstate(jj) == ACTIVE
    if length(hdp.dp{jj}.beta)~=hdp.base.numclass+1
      error(sprintf('length(hdp.dp{%d}.beta)!=hdp.base.numclass+1',jj));
    end
    if any(hdp.dp{jj}.beta < 0) 
      error(sprintf('hdp.dp{%d}.beta<0',jj));
    end
    if abs(sum(hdp.dp{jj}.beta) - 1) > 1e-8
      error(sprintf('sum(hdp.dp{%d}.beta)!=1',jj));
    end
  end
end

% Check classqq
for cc = 1:hdp.base.numclass+1
  ni = qq_numitems(hdp,cc);
  if ni~=0
    error(sprintf('hdp.base.classqq(:,%d) inconsistent with data',cc));
  end
end

% Check classnd, classnt, datacc
for jj = hdp.numdp:-1:1
  if hdp.dpstate(jj) == ACTIVE | hdp.dpstate(jj) == FROZEN
    oo = ones(size(hdp.dp{jj}.datacc));
    nd = full(sparse(oo,hdp.dp{jj}.datacc,oo,1,hdp.base.numclass+1));
    hdp.dp{jj}.classnd = hdp.dp{jj}.classnd - nd;

    if any(hdp.dp{jj}.classnd ~= 0)
      error(sprintf('hdp.dp{%d}.classnd inconsistent with datacc',jj));
    end
    pp = hdp.ppindex(jj);
    if pp == 0
      ii = find(hdp.dp{jj}.classnd>0);
      if any(hdp.dp{jj}.classnt(ii)~=1)
        error(sprintf('hdp.dp{%d}.classnt!=1',jj));
      end
      if hdp.dp{jj}.classnt(hdp.base.numclass+1)~=0
        error(sprintf('hdp.dp{%d}.classnt(hdp.base.numclass+1)!=0',jj));
      end
    else
      hdp.dp{pp}.classnd = hdp.dp{pp}.classnd - hdp.dp{jj}.classnt;
    end
  end
end
