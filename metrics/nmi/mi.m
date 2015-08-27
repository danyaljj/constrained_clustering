function score = mi(cmat)

% function score = mi(cmat)
%   compute normalized mutual information (NMI) from the confusion matrix
%   'cmat' and return the value in 'score'

nh = sum(cmat,1);
nl = sum(cmat,2);
[s1,s2] = size(cmat);
nlh = repmat(nl,[1,s2]) .* repmat(nh,[s1,1]); %nlh = nl * nh;
ind = find(nlh > 0);

nlh(ind) = cmat(ind) ./ nlh(ind);
n = sum(nl);
nlh = nlh * n;

ind = find(nlh>0);
nlh(ind) = log(nlh(ind));
nlh = cmat .* nlh;
score = sum(nlh(:));

ind = find(nh>0);
nh(ind) = nh(ind) .* log(nh(ind)/n);
ind = find(nl>0);
nl(ind) = nl(ind) .* log(nl(ind)/n);
tmp = sqrt(sum(nh)*sum(nl));
if tmp == 0
    score = 0;
else
    score = score / tmp;
end

return;
