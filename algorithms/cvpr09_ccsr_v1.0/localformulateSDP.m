function [F0, FI, c] = localformulateSDP(S, D, b)
% formulate SDP problem
% each FI that corresponds to the LMI for the quadratic cost function has
% precisely 2*D^2 nonzero elements. But we need only D^2 storage for
% indexing these elements since the FI are symmetric
tempFidx = zeros(D^2, 3);
dimF = (D^2+1) + D;
idx= 0;
for col=1:D
    for row=col:D
        idx = idx+1;
        lindx1 = sub2ind([D D], row, col);
        lindx2 = sub2ind([D D], col, row);
        tempFidx(:,1) = [1:D^2]';
        tempFidx(:,2) = D^2+1;
        if col==row
            tempFidx(:,3) = S(:, lindx1) ;
            FI{idx} = sparse([tempFidx(:,1); ...  % for cost function
                                tempFidx(:,2); ... % symmetric
                                row+D^2+1 ... % for P being p.s.d
                              
                            ], ...
                            [tempFidx(:,2); ...  % for cost function
                                tempFidx(:,1); ... % symmetric
                                row+D^2+1; ... % for P being p.s.d
                               
                            ],...
                            [tempFidx(:,3); ... % for cost function
                                tempFidx(:,3); ... % symmetric
                                1;                  % for P being p.s.d
                               
                            ], dimF, dimF);
        else
            
            tempFidx(:,3) = S(:, lindx1) + S(:, lindx2);
            FI{idx} = sparse([tempFidx(:,1); ...  % for cost function
                                tempFidx(:,2); ... % symmetric
                                row+D^2+1; ... % for P being p.s.d
                                col+D^2+1; ... % symmetric
                            ], ...
                            [tempFidx(:,2); ...  % for cost function
                                tempFidx(:,1); ... % symmetric
                                col+D^2+1; ... % for P being p.s.d
                                row+D^2+1; ... % being symmetric
                            ],...
                            [tempFidx(:,3); ... % for cost function
                                tempFidx(:,3); ... % symmetric
                                1;                  % for P being p.s.d
                                1;                  % symmetric
                            ], dimF, dimF);
            
        end
    end
end
idx=idx+1;
% for the F matrix corresponding to t
FI{idx} = sparse(D^2+1, D^2+1, 1, dimF, dimF);

% now for F0
F0 = sparse( [[1:D^2]], [[1:D^2]], [ones(1, D^2)], dimF, dimF);

% now for c
b = reshape(-b, D, D);
b = b*2 - diag(diag(b)); 
c = zeros(idx-1,1);
kdx=0;
%keyboard;
for col=1:D
    for row=col:D
      kdx = kdx+1;
      c(kdx) = b(row, col);
    end
end
%keyboard;
c = [c; 1]; % remember: we use only half of P
return;
