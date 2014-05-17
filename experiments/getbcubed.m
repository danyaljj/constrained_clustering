function [precision , recall ] = getbcubed(assignmentGold, assignmentPred) 

numPred = 0; 
numGold = 0; 
numBoth = 0; 
for i = 1:length(assignmentGold)
    
    %     assignmentGold
%     length(assignmentGold)
    for j = 1:length(assignmentGold)
        if( i == j ) 
            continue; 
        end
        if( assignmentGold(i) == assignmentGold(j) && assignmentPred(i) == assignmentPred(j) ) 
            numBoth = numBoth + 1; 
        end
        if( assignmentGold(i) == assignmentGold(j)  ) 
            numGold = numGold + 1; 
        end
        if(  assignmentPred(i) == assignmentPred(j) ) 
            numPred = numPred + 1; 
        end        
    end
end

% numBoth 
% numGold 
% numPred 

precision = numBoth / numGold;
recall = numBoth / numPred; 

