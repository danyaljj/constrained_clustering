function ww = randdir(aa,normdim);

ww = randgamma(aa);

if length(aa) <= 1
  ww = 1;
  return;
end

if nargin == 1
  normdim = find(size(aa)>1);
  normdim = normdim(1);
end

ndim = ndims(aa);
index = struct('type','()','subs',{repmat({':'},1,ndim)});
index.subs{normdim} = ones(1,size(aa,normdim));

ww = ww./subsref(sum(ww,normdim),index);
