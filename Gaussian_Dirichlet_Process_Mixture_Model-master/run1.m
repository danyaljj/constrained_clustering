clear all; close all; clc; 

% sampled data:
mean1 = [2,-1];
cov1 = [1,0.1; 0.1,1];
mean2 = [8,3];
cov2 = [1 .2; 0.2,1];
X = [mvnrnd(mean1, cov1, 200); mvnrnd(mean2, cov2, 200)];

T = 50;  % maximum number of clusters
[gamma, phi, m, beta, s, p] = variational_dpm(X, 20, T, 1);
[maxVal, clusters] = max(phi);
for t = 1:T 
    xt = X(clusters == t, :);
    if size(xt) ~= 0
        disp( ['T = ' num2str(t) ' size(xt,1) = '  num2str(size(xt,1))  '  m(t,:) ' num2str(m(t,:)) ])
    end 
end

plot(X(1:200, 1), X(1:200, 2), 'xr')
hold on; 
plot(X(201:end, 1), X(201:end, 2), 'xb')

% % sample output : 
% T = 4 size(xt,1) = 221  m(t,:) 7.5048      2.7682
% T = 50 size(xt,1) = 179  m(t,:) 2.0489     -1.1736
plot(7.5048,2.7682,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
plot(2.0489,-1.1736,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
