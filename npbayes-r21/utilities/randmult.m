function ii = randmult(pp,normdim);

% pp need not be normalized

ss    = size(pp);
if nargin == 1
  normdim = find(ss>1);
  if isempty(normdim)
    ii = 1; 
    return
  elseif length(normdim)==1
    pp = cumsum(pp);
    ii = 1+sum(rand*pp(end)>pp);
    return
  end
  normdim = normdim(1);
end

ndim                = ndims(pp);
col                 = {':'};
index               = struct('type','()','subs',{col(1,ones(1,ndim))});

pp                  = cumsum(pp,normdim);

tt                  = ss;
tt(normdim)         = 1;
index.subs{normdim} = ss(normdim);
rr                  = rand(tt) .* subsref(pp,index);

index.subs{normdim} = ones(1,ss(normdim));
ii                  = 1+sum(subsref(rr,index)>pp,normdim);

