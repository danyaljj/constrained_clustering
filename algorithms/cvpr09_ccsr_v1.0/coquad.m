function [B, b] = coquad(Q,M,C)

% Formulate eq.(18) in paper
% Li,et al., CVPR 2009, "Constrained Clustering via Spectral Regularization" 
% Please report problems to Zhenguo Li at zgli@ee.columbia.edu

[Npts m] = size(Q);
B = zeros(m^2);
b = zeros(m^2,1);

for i = 1:Npts
    U = Q(i,:)'*Q(i,:);
    s = U(:);
    B = B + s*s';
    b = b + s;
end

for k = 1 : size(M,1)
    i = M(k,1);
    j = M(k,2);
    U = Q(j,:)'*Q(i,:);
    s = U(:);
    B = B + s*s';
    b = b + s;
end

for k = 1 : size(C,1)
    i = C(k,1);
    j = C(k,2);
    U = Q(j,:)'*Q(i,:);
    s = U(:);
    B = B + s*s';
end

b = -2 * b;