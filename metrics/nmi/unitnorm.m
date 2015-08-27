function dat = unitnorm(dat,dim)

% dat = unitnorm(dat, dim)
%   normalize dat along given dimension (default dim=1, i.e. each
%   column gets normalized to unit length in L2-norm)

if nargin < 2, dim = 1; end

nd = sum(dat.^2, dim);
I = find(nd>0);
nd(I) = 1 ./ sqrt(nd(I));
nd = diag(sparse(nd));
if dim == 2
    dat = nd * dat;
else
    dat = dat * nd;
end

return;

