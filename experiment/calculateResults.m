function [results] = calculateResults(dataAll, assignmentAll, titles, trueLabels)

results{1,1} = 'Method';
results{1,2} = 'P';
results{1,3} = 'R';
results{1,4} = 'F';
results{1,5} = 'AdjRnd';
results{1,6} = 'NMI';
for i = 1:length(dataAll)
    disp([ 'i = ' num2str(i) ' - name = ' titles{i} ])
    if( ~isempty( assignmentAll{i} )  )
        [vnmi, p, r, f1, adjrnd] = evaluateAll( assignmentAll{i}, trueLabels);
    else
        vnmi = 0; p = 0; r = 0; f1 = 0; adjrnd = 0; 
    end 
    results{i+1,1} = titles{i};
    results{i+1,2} = p;
    results{i+1,3} = r;
    results{i+1,4} = f1;
    results{i+1,5} = adjrnd;
    results{i+1,6} = vnmi;
end
end