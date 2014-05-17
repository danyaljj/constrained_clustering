function mm = logmeanexp(xx, normdim);

if nargin == 1
  normdim = find(size(xx)>1);
  normdim = normdim(1);
end

ndim = ndims(xx);
index = struct('type','()','subs',{repmat({':'},1,ndim)});
index.subs{normdim} = ones(1,size(xx,normdim));

xmax = max(xx,[],normdim);

mm = xmax + log(mean(exp(xx-subsref(xmax,index)),normdim));

