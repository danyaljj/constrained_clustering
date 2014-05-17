function check_ttindex(cpindex,ttindex);

numconparam = max(cpindex);

for ii = 1:numconparam
  jj = find(cpindex==ii);
  tt = ttindex(jj);
  if any(sort(tt)~=1:length(tt))
    error('Not valid ttindex');
  end
end

