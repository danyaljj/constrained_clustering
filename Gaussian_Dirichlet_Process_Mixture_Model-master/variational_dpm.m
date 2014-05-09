function [gamma, phi, m, beta, s, p] = variational_dpm(x, iter, T, alpha)

if isempty(alpha)
    alpha = 0.1;
end
if isempty(T)
    T = 30;
end
if isempty(iter)
    iter = 30;
end

gamma = ones(T-1, 2);

phi = rand( [T size(x,1)] );
phi = phi ./ repmat( sum(phi), [size(phi,1), 1] );

s1 = 1;
s2 = 1;

m0 = zeros(1, size(x,2));
beta0 = 1;
p0 = size(x,2) + 2;
s0 = eye(size(x,2));

m = rand( [T size(x,2)] );
beta = ones([T 1]);
s = rand( [T size(x,2) size(x,2)] );
p = 2*ones( [T 1] );

for i = 1:iter
    gamma(:, 1) = 1 + sum(phi(1:T-1, :), 2);
    tmp = cumsum(phi(end:-1:2, :));
    phi_cum = tmp(end:-1:1, :);
    gamma(:,2) =  alpha + sum(phi_cum, 2);
    
    [maxVal, clusters] = max(phi);
    
    for t = 1:T
        xt = x(clusters == t, :);
        if size(xt, 1) == 0
            m(t, :) = m0;
            beta(t) = beta0;
            p(t) = p0;
            s(t,:,:) = s0;
            continue
        end
        meanxt = mean(xt);
        m(t, :) = 1. * (size(xt,1) * meanxt + beta0 * m0) / (size(xt,1) + beta0);
        beta(t) = beta0 + size(xt,1);
        p(t) = p0 + size(xt,1);
        s(t, :, :) = s0 +  (xt - repmat(meanxt, [size(xt,1) 1]))' *  (xt - repmat(meanxt, [size(xt,1) 1])) / size(xt,1) + ...
            (size(xt,1) * beta0 / (size(xt,1) + beta0) ...
            * ((meanxt - m0)' * (meanxt - m0)));
    end
    
    sumGamma = gamma(:, 1) + gamma(:, 2);

    Eq_logv1 = digamma(gamma(:, 1)) - digamma(sumGamma);

    Eq_logv2 = digamma(gamma(:, 2)) - digamma(sumGamma);
    
    cumsum_Eq_logv2 = cumsum(Eq_logv2);
    
    Eq_pi = zeros(1,T);
    Eq_pi(1) = exp(Eq_logv1(2));
    Eq_pi(2:end-1) = exp(Eq_logv1(2:end) + cumsum_Eq_logv2(1:end-1));
    Eq_pi(end) = 1 - sum(Eq_pi(1:end));
    
    Eq_mu = m;
    Eq_sigma = s; 
    for t = 1:T
        tmp = zeros([size(x,1) 1]);
        for ii = 1:size(x,1)
            tmp(ii) = gaussian_log_pdf(x(ii,:), Eq_mu(t,:), reshape(Eq_sigma(t,:,:), [size(Eq_mu(t,:),2) size(Eq_mu(t,:),2)]));
        end
        phi(t, :) = tmp;
    end
    phi = exp(phi - repmat(max(phi), [ size(phi,1) 1]));
    phi = (phi .* repmat(reshape(Eq_pi, T, 1),[1, size(phi,2)]) );
    phi = phi ./ repmat( sum(phi), [size(phi,1) 1] );	%normalize
end