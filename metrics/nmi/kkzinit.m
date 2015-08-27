function mu = kkzinit(dat, kc)

% mu = kkz(dat, kc)
% initialize according to Katsavounidis, Kuo, and Zhang (1994)

mu(:,1) = unitnorm(sum(dat,2),2);
[y, i] = min(mu'*dat); % take the first centroid to be the doc
                       % most distant to global mean
%[y, i] = max(sum(dat>0,1));
mu(:,1) = dat(:,i);

[d, n] = size(dat);
dist = zeros(1,n);
ind = zeros(1,n);

k = 1;
while k < kc
    dd = mu(:,k)' * dat;
    [dist, I] = max([dd; dist],[],1);
    ind(I==1) = k;
    [y, i] = min(dist);
    k = k + 1;
    mu(:,k) = dat(:,i);
end
    
return;
