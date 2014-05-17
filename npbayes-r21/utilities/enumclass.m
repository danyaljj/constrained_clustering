function numclass = enumclass(alpha,numdata);

numclass = alpha*sum(1./(alpha-1+(1:numdata)));
