function check_ppindex(numdp,ppindex);

if any(ppindex<0)|any(ppindex>=1:numdp)|any(ppindex~=ceil(ppindex))
  error('Not valid ppindex');
end

