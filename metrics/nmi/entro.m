function e = entro(x)

% e = entro(x)
%
% calculate the entropy of a non-negative vector x, which will be
% normalized first so that sum(x) = 1

if min(size(x)) ~= 1
	error('input x must be a vector');
end
if ~isempty(find(x < 0))
	error('input x must be non-negative');
end

p = x / sum(x);
p(p==0) = 1;
lp = log2(p);

if size(p,1) == 1
    e = - p * lp';
else
    e = - lp' * p;
end

return;

