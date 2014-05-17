function qq = Multinomial_additems(hh,qq,cc,ss);

sc = sparse(ss,cc,ones(size(ss)),size(qq,1),size(qq,2));
qq = qq + sc;

