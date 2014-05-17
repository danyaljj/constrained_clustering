function beta = randstick(alpha,numclass);

one = ones(1,numclass-1);
zz = [randbeta(one, alpha*one) 1];
beta = zz .* cumprod([1 1-zz(1:numclass-1)]);
