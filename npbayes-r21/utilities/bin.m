function bb = bin(cc,loc);

ss = size(loc);
mm = max(loc);
cc = cc(:);
oo = ones(1,length(cc));
bb = sparse(oo,cc,oo,1,mm);
bb = reshape(bb(loc),ss);
