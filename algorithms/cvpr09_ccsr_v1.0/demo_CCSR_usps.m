% Codes for the paper
% "Constrained Clustering via Spectral Regularization", CVPR 2009
% Zhenguo Li, Jianzhuang Liu, and Xiaoou Tang.
% Written by Zhenguo Li, zgli@ee.columbia.edu
% Version 1.0, Dec. 01, 2010

clc,clear,close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load data

load usps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. parameters

k = length(unique(labels)); % number of clusters
Npts = size(X,1); % number of points
dm = squareform(pdist(X)); % distance matrix
r = averagekmin_dm(dm,20);
r1 = linspace(0.1*r, r, 5);
r2 = linspace(r, 10*r, 5);
sigma = unique([r1,r2]); % set of scale factors in graph construction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. generate pairwise constraints

[ML CL] = genPWC(labels,20,20); % generate pairwise constraints randomly

err_sigma = zeros(length(sigma),1);
for i_sigma = 1:length(sigma)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. compute normalized graph Laplacian
    
    L = spnlaplacian_dm(dm,sigma(i_sigma),20); % sparse normalized Laplacian
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. compute the first m eigenvectors of normalized graph Laplacian
    
    m = 15; % number of eigenvectors used in the paper
    opts.disp = 0;
    [Q,E] = eigs(L,m,'sm',opts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6. formuate the convex quadratic semidefinite program
    [A, b] = coquad(Q,ML,CL);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 7. formuate and solve the semidefinite program
    
    S = sqrtm(A); S = real(S); S = (S + S')/2; % matrix squared root
    
    % symmetrize b, necessary due to formulateSDP below
    b = reshape(b,[m,m]);b = (b+b')/2;
    b = b(:);
    
    [AA, bb, cc] = formulateSDP(S, m, -b); % formulate the SDP
    K.s = m^2 + 1 + m;
    [xx, yy, zz, info] = csdp(AA, bb, cc,K); % solve the SDP
    yy = -yy; % the negative of yy is our solution
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 8. obtain the low-dimensional embedding and call kmeans
    
    M = getY(yy,m);
    P = sqrtm(M); P = real(P); P = P + P';
    Y = Q * P;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9. call kmeans and show results
    
    res = kmeans(Y,k,'Replicate',10);
    
    err_sigma(i_sigma) = get_error_rate(labels,res);
    
    sprintf('%d out of %d, the error rate is: %f\n', i_sigma, length(i_sigma), err_sigma(i_sigma)),    
end

sprintf('the best error rate: %f\n',min(err_sigma)),