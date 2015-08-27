function mu = perturbinit(dat, kc)

d = size(dat,1);
mu = repmat(sum(dat,2),[1,kc]);
mu = mu .* (1+(rand(d,kc)-0.5)/10);
mu = unitnorm(mu,1);

return;

