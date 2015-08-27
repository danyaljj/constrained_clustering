function partition = bpart(ll,m)

% partition = bpart(ll,m)
%
% Balanced partitioning of a log-likelihood matrix (ll, size n by kc).
% The second argument (m) is optional and, if given, means that
% at least m samples must be in each cluster. m=0 corresponds to
% a completely balanced partition.

if nargin < 2, m = 0; end

[n,k] = size(ll);

% balanced re-partition
pn = 1:n;
partition = zeros(n,1);
noti = 1:k;

if m > 0
    % guarantee each cluster to have at least m samples
    for i = 1 : k
        noti = setdiff(noti,i);
        if i < k
            dll = ll(pn,noti) - repmat(ll(pn,i),[1,length(noti)]);
        else
            dll = - ll(pn,i);
        end
        [oll, order] = sort(max(dll,[],2));
        tmp2 = pn(order(1:m));
        partition(tmp2) = i;
        pn = setdiff(pn, tmp2);
    end
    % use ML assignment for rest samples
    [oll, order] = max(ll(pn,:),[],2);
    partition(pn) = order;
else
    for i = 1 : k
        %disp(sprintf('i = %d',i));
        nci = round(i*n/k) - round((i-1)*n/k);
        noti = setdiff(noti,i);
        if i < k
            dll = ll(pn,noti) - repmat(ll(pn,i),[1,length(noti)]);
            [oll, order] = sort(max(dll,[],2));
            tmp2 = pn(order(1:nci));
        else
            tmp2 = pn;
        end
        partition(tmp2) = i;
        pn = setdiff(pn, tmp2);
    end
end

return;

