function [h] = plotExperiments(dataAll, assignmentsAll, titles, height, width)
h = figure;
hold on;
colormap([1 0 .5;   % magenta
           0 0 .8;   % blue
           0 .6 0;   % dark green
           .3 1 0]); % bright green
dotsize = 12;
DataSize = length(dataAll); 
for i=1:DataSize
    subplot(height,width,i);
    data = dataAll{i}; 
    assignments = assignmentsAll{i}; 
    scatter(data(:,1), data(:,2), dotsize, assignments); axis equal;
    t = titles{i};
    title(t);
end 
hold off;
end