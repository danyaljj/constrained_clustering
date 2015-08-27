function out = cm(class, cluster)

% function cmat = cm(class, cluster)
%   Compute the confusion matrix (cmat) given class labels (class)
%   and cluster IDs (cluster)

if size(class,1)==1
    cla = class';
else
    cla = class;
end
n = length(class);
a = spconvert([(1:n)' cla ones(n,1)]);

if min(size(cluster)) == 1
    if size(cluster,1)==1
        clu = cluster';
    else
        clu = cluster;
    end
    b = spconvert([(1:n)' clu ones(n,1)]);
    b = b';
else
    if size(cluster,1) == n
        b = cluster';
    else
        b = cluster;
    end
end

out = full(b * a);

return;
