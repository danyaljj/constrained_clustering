function check_dpindex(numdp,dpindex);

if any(dpindex<=0) | any(dpindex>numdp) | any(dpindex~=ceil(dpindex)) | ...
   length(dpindex)~=length(unique(dpindex))
  error('Not valid dpindex');
end

