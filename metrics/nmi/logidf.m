function dat = logidf(dtm)

% DAT = LOGIDF(DTM)
%
% To weight each word using log of its inverse document frequency (IDF)
% this is usually used for vMF models. DTM is the input document-term
% matrix and DAT is the transformed data matrix. Both have a size of n by d.

n = size(dtm,1);
idf = sum(dtm>0,1);
I = find(idf>0);
idf(I) = n ./ idf(I);
idf(I) = sparse(log(idf(I)));
sidf = diag(idf);
dat = dtm * sidf;

return;

