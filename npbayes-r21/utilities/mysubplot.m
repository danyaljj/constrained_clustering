function mysubplot(i,j,k,m);

if i>0 & j>0,

 if nargin==3, m=0; end;

 if length(k)==1,
  jk = rem(k-1,j);
  ik = i-1-fix((k-1)/j);
 elseif length(k)==2,
  jk = k(2)-1;
  ik = i-(k(1));
 end;

 subplot('position',[jk/j+m .95*ik/i+m 1/j-2*m .95/i-2*m]);

else
 error('At least one subplot');
end;
