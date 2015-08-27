clear;

number = '5'; 
fileNamesAll = {'Two-Spirals', ...  
                'Cluster-In-Cluster',  ... 
                'crescentfullmoon', ... 
                'Gaussian-Mixtures', ... 
                'Half-Kernel', ... 
                'outlier'};
height = 6;
width = 8;
dotsize = 12;
h = figure;
colormap([1 0 .5;   % magenta
    0 0 .8;   % blue
    0 .6 0;   % dark green
    .3 1 0]); % bright green
hold on;
last_i = 0;     
for fileName = fileNamesAll
    fileName
    load([fileName{1} number])
    %h = plotExperiments(dataAll, assignmentAll, titles, 4, 5, title);
    % title(t)
    assignmentAll(5) = [];
    titles(5) = [];
    dataAll(5) = [];
    DataSize = length(dataAll); 
    for i=1:DataSize
        subplot(height,width,i + last_i);
        data = dataAll{i};
        assignments = assignmentAll{i};
        scatter(data(:,1), data(:,2), dotsize, assignments); axis equal;
        t = titles{i};
        title(t)
    end
    last_i = last_i + DataSize; 
end
