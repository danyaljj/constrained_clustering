function numtable = randnumtable(weights,maxtable);

numtable = zeros(size(maxtable));
[B I J]  = unique(maxtable(:));

weights = log(weights);

for ii = 1:length(B)
  maxtable = B(ii);
  if maxtable > 0
    mm = 1:maxtable;
    stirnum = stirling(maxtable);
    for jj = find(J==ii)'
      clik = mm .* weights(jj);
      clik = cumsum(stirnum .* exp(clik-max(clik)));
      numtable(jj) = 1+sum(rand*clik(maxtable) > clik);
    end
  end
end


