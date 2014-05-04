function [Phi, C, W] = gaussian_EM(X, n_comp, n_iter)
% USAGE:
% G4 = [-12, 1.0;
%       -4, 1.0;
%        0, 1.0
%        9, 1.0];
% W4 = [0.2 0.4 0.2 0.2];
% n=50;
% X=zeros(n,1);
% for i =1:n
%     g = G4(mnrnd(1,W4)==1,:);
%     X(i) = normrnd(g(1), g(2));
% end
% [Phi, C, W] = gaussian_EM(X, 4, 10);

% Proper version:
n = length(X);

% randomized normalized assignment vector
% C = rand(n, n_comp);
% sums = sum(C,2);
% cs = sums(:,ones(1,n_comp));
% C = C./cs;

% kmeans inital assignment vector
IDX = kmeans(X,n_comp);
C=zeros(length(X), n_comp);
for i=1:n_comp
    C(IDX==i,i) = 1;
end

% Param vector
Phi = zeros(n_comp, 2);
sc = zeros(n_comp, 1); % for summing columns of C
for c=1:n_comp
    sc(c) = sum(C(:,c));
    Phi(c,1) = sum(C(:,c) .* X) / sc(c);
    Phi(c,2) = sum(C(:,c) .* (X - Phi(c,1)).^2) / sc(c);
    % hack to protect mean and variance
    if abs(Phi(c,2)) < 1e-75 || isnan(Phi(c,2)); Phi(c,2) = 0.1; end
end

% Weight vector
W = ones(n_comp, 1)/n_comp;

for s=1:n_iter
    % E step: evaluate component assignment probabilities given current
    % mus and sigmas
    for i=1:length(X)
        ls = zeros(n_comp, 1); % likelihoods
        for c=1:n_comp
            ls(c) = W(c) * normpdf(X(i), Phi(c,1), Phi(c,2));
        end
        C(i,:) = ls/sum(ls); % normalize
        if sum(isnan(C(i,:)))>0
            error('isnan in c');
        end
    end
    
    % M step: compute ML mean and variance given current assignments
    sc = zeros(n_comp, 1); % for summing columns of C
    for c=1:n_comp
        sc(c) = sum(C(:,c)); if sc(c)==0; sc(c)=1; end
        Phi(c,1) = sum(C(:,c) .* X) / sc(c);
        Phi(c,2) = sum(C(:,c) .* (X - Phi(c,1)).^2) / sc(c);
        
        % hack to protect mean and variance
        if abs(Phi(c,2)) < 1e-75 || isnan(Phi(c,2)); Phi(c,2) = 0.1; end
    end
    if sum(isnan(Phi))>0
        error('isnan phi');
    end
    W = sc/sum(sc);
end

% plot results
clf;
hold on;
hist(X,100);
xrange=-15:0.1:15;
y = zeros(1,length(xrange));
for i=1:n_comp
    y = y + 20 * W(i) * normpdf(xrange, Phi(i,1), Phi(i,2));
end
plot(xrange, y);

end