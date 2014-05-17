function alpha = randconparam(alpha,numdata,numclass,aa,bb,numiter);

% Escobar and West's method for single Gamma.

if nargin == 5
  numiter = 1;
end

for ii = 1:numiter
  xx = randbeta(alpha+1,numdata);

  weights = [(aa+numclass-1)/(bb-log(xx)) numdata];
  zz = randmult(weights/sum(weights));

  if zz==1
    alpha = randgamma(aa+numclass,1./(bb-log(xx)));
  else
    alpha = randgamma(aa+numclass-1,1./(bb-log(xx)));
  end
end

