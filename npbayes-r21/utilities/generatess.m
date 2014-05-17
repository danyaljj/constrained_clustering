function ss = generatess(data);

numdoc = size(data,1);
ss = cell(1,numdoc);

for dd = 1:numdoc
  word = find(data(dd,:));
  freq = data(dd,word);
  
  cs = [0 cumsum(freq)];
  ss{dd} = zeros(1,sum(data(dd,:)));
  
  for ii=1:length(word)
    ss{dd}(cs(ii)+1:cs(ii+1)) = word(ii);
  end

end
