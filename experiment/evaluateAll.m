function [vnmi, p, r, f1, ajrnd] = evaluateAll( assignment, trueLabels )
    min1 = min(min(assignment)); 
    if min1 <= 0
        assignment = assignment + 1 - min1; 
    end
    min1 = min(min(trueLabels));    
    if min1 <= 0  
        trueLabels = trueLabels + 1 - min1; 
    end
    [p, r] = getbcubed(trueLabels, assignment); 
    %vnmi = nmi(trueLabels, assignment); 
    cmatrix = cm(trueLabels,assignment);
    vnmi = mi(cmatrix);
    f1 = 2*p*r / (p+r);
    ajrnd = adjrand(assignment, trueLabels); 
end 