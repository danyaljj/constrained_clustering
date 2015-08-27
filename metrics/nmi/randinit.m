function mu = randinit(dat, kc)

n = size(dat,2);
r = mod(randperm(n),kc)+1;
rp = spconvert([(1:n)' r' ones(n,1)]);
mu = dat * rp;
mu = unitnorm(mu,1);

return;

