function alpha = randconparam2(alpha,numdata,numclass,aa,bb);

% less stable for small alpha's as xx gets very close to zero.

xx = randbeta(alpha,numdata);

alpha = randgamma(aa+numclass,1./(bb-log(xx)));

