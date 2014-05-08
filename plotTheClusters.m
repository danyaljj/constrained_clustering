function plotTheClusters(X, clusterSize,  pointsInCluster)

% Xtmp = X(Y ==1, :);
plot(X(:, 1), X(:, 2), 'xr')
hold on;
% Xtmp = X(Y ==-1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')
for i = 1:clusterSize
    plot(X(pointsInCluster(i),1), X(pointsInCluster(i),2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
end

end