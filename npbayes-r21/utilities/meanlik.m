function lik = meanlik(lik);

lik = logmeanexp(-logmeanexp(-lik,3),2);
