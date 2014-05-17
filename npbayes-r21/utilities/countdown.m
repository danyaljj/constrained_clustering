function countdown(fid,preamble,numiter,totiter);

persistent starttime prevtime totaliter curiter

if nargin==4
  starttime = toc;
  prevtime  = toc;
  totaliter = totiter;
  curiter   = numiter;
  fprintf(fid,preamble,totiter);
else
  curiter   = curiter + numiter;
end

if toc-prevtime > 5
  tt = toc-starttime;
  ss = sprintf('time %1.1f ETC %1.1f',tt,tt/curiter*totaliter);
  fprintf(fid,[preamble ss]);
  prevtime = toc;
end

