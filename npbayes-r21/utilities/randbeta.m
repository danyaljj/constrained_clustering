function xx = randbeta(aa,bb);

shape = size(aa);

xx = randgamma(cat(2,aa(:),bb(:)));
xx = reshape(xx(:,1)./sum(xx,2),shape);

