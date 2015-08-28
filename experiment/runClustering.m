function [] = runClustering(X, Y, k, experimentName, rate, p, verboseOutput)
% X: The input data 
% Y: The output labels 
% k: the desired number of clusters 
% experimentName: the name of the experiement; for logging results. 
% rate: the size of the pairwise constraints over the total possible 
%       number of constraints. 
% p: a value in [0, 1]. Thi is the confidence on quality of the constraints. 
% When p = 1, the constraints are high quality (noise-free). As `p` gets 
% closer to zero, it is more likely that the constraints will be fliped 
% (hence noisier). 
% verboseOutput: a boolean input. If this is true, the code will visuzlize 
% the output and will save it and the results on disk. Note that the 
% visualization makes sense only if the data is 2D 

iter = 1;
name1 = [experimentName '_iter=' num2str(iter)];
name1 = strrep(name1, '.', '_');

% checking if there is a result on disk with the same experiment name. 
% if there is one, we open it, and add our results into it, and save it
% back. 
if ~exist([name1 '.mat'])
    disp(['skipping the experiment '  experimentName])
    dataAll = {};
    assignmentAll = {};
    titles = {};
    calculationTime = {};
        
    % make E matrix: the matrix containing pairwise constraints 
    % since different codes have different input standards, here I am
    % creating multiple encodings of the constraint matrix. 
    E = zeros(size(Y, 1), size(Y, 1));
    Checked = zeros(size(Y, 1), size(Y, 1));
    randSize = rate * size(Y, 1) * size(Y, 1);
    iterAll = 1;
    ML = [];
    CL = [];
    C_m = [];
    while(1)
        i1 = randi(size(Y, 1));
        i2 = randi(size(Y, 1));
        if i1 == i2
            continue;
        end
        if Checked(i1, i2) == 0
            correct = (rand() < p);
            Checked(i1, i2) = 1;
            if Y(i1) == Y(i2)
                if correct == 1
                    E(i1, i2) = 1;
                    E(i2, i1) = 1;
                    ML = [ML; i1 i2]; % TODO: should we add the reverse too?
                    ML = [ML; i2 i1];
                    C_m = [C_m; i1 i2 1];
                    C_m = [C_m; i2 i1 1];
                else
                    E(i1, i2) = -1;
                    E(i2, i1) = -1;
                    CL = [CL; i1 i2];
                    CL = [CL; i2 i1];
                    C_m = [C_m; i1 i2 -1];
                    C_m = [C_m; i2 i1 -1];
                end
            else
                if correct == 1
                    E(i1, i2) = -1;
                    E(i2, i1) = -1;
                    CL = [CL; i1 i2];
                    CL = [CL; i2 i1];
                    C_m = [C_m; i1 i2 -1];
                    C_m = [C_m; i2 i1 -1];
                else
                    E(i1, i2) = 1;
                    E(i2, i1) = 1;
                    ML = [ML; i1 i2];
                    ML = [ML; i2 i1];
                    C_m = [C_m; i1 i2 1];
                    C_m = [C_m; i2 i1 1];
                end
            end
            iterAll = iterAll + 1;
        end
        if( iterAll > randSize)
            break;
        end
    end
    matConst = E + eye(size(Y, 1));
    
    T = mean(X);
    [dist, ind] = sort( sqrt(sum((repmat(T,size(X,1),1)-X).^2,2)), 'descend' );
    lambda = dist(k);
    E2 = C_m;
    E2(:,3) = (E2(:,3) + 1)/2;    
else
    disp('loading the file . . . ')
    load([name1 '.mat'])
end

%% true labels
name = 'True labels';
if ~containsMethod(name, titles)
    ind = length(dataAll);
    dataAll{ind+1} = X;
    assignmentAll{ind+1} = Y;
    titles{ind+1} = name;
end

%% k-means
name = 'K-means';
if ~containsMethod(name, titles)
    disp('K-means ....');
    tic;
    [~, ~, assignment]= kmeans2(X, k);
    calculationTime{ind+1} = toc;
    ind = length(dataAll);
    dataAll{ind+1} = X;
    assignmentAll{ind+1} = assignment;
    titles{ind+1} = name;
    disp('kmeans is done!!')
end

%% rdp-means
if 1 % easy way to activate/diactivate this part of code 
    name = 'RDP-means';
    if ~containsMethod(name, titles)
        tic;
        try
            [~, ~, assignment_out, ~, ~, ~, centroindsAll] = rdpmeans(X, lambda, E, 2, .001);
        catch
            assignment_out = [];
            disp('Error handled')
        end
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = assignment_out;
        titles{ind+1} = name;
    end
end


%% LCVQE
if 1
    name = 'LCVQE';
    if ~containsMethod(name, titles)
        tic;
        try
            [idx, ~, ~, ~, ~] = lcvqe(X, k, C_m);
        catch 
            idx = [];
            disp('Error handled')
        end
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = idx;
        titles{ind+1} = name;
    end
end

%% dp-means:
if 1
    name = 'DP-means';
    if ~containsMethod(name, titles)
        tic;
        disp('DP-means ....');
        disp(['experiment name  = ' experimentName])
        [~, ~, assignment, clusterSize]= dpmeans(X, lambda);
        ind = length(dataAll);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
        % save('data')
    end
end

%% constrained k-means
if 0 % this not active since this algorithm is does not always return 
     % responses, especially when data is noisy. 
    [centroid, pointsInCluster, assignment] = constrained_kmeans(X, E, k);
    ind = length(dataAll);
    dataAll{ind+1} = X;
    if length(assignment) < size(X,1)
        maxComp = max(assignment);
        diffSize = size(X,1) - length(assignment);
        assignment = [assignment; randint(diffSize, 1, maxComp) + 1];
    end
    assignmentAll{ind+1} = assignment;
    titles{ind+1} = 'Constrained K-means';
end 

%% MPCKMeans
% Based on a java code by Bilenko et al. (2004). The following function 
% call is to a matlab wrapper of java function. 
if 1
    name = 'MPCKMeans';
    if ~containsMethod(name, titles)
        tic;
        [assign1, cent1] = runMPCKMeans(X, Y, k, C_m);
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = double(assign1);
        titles{ind+1} = name;
    end
end

%% Co1SC, 2012: constrained 1-spectral clustering 
if 1
    name = 'Co1SC-NC';
    if ~containsMethod(name, titles)
        disp('Co1SC, 2012: constrained 1-spectral ...');
        % W:      n x n similarity matrix, where n is the number of data points.
        %           It has to be a sparse symmetric matrix, with zeros on the diagonal.
        % deg:    n x 1 degree vector (can also be seen as vertex weights)
        %           for Normalized cut it is sum(W,2)
        %           for Ratio cut, it is ones(size(W,1),1)
        % k:      number of clusters
        % ML:     m x 2 martix, specifying m must-link pairs
        % CL:     p x 2 matrix, specifying p cannot-link pairs
        tic;
        try
            W_nonsparse = squareform(pdist(X));
            W = sparse(W_nonsparse);
            vertex_weights = sum(W,2); % this choice corresponds to normalized cut
            [~, assignment, viols] = cosc(W, vertex_weights, k, ML, CL);
            % Faster, but approximate form: 
            % [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000);
        catch
            assignment = [];
            disp('Error handled')
        end
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
    end
end

% 1-clustering, with different definition of Laplacian. 
if 1
    name = 'Co1SC-RC';
    if ~containsMethod(name, titles)
        tic;
        try
            vertex_weights = ones(size(W,1),1); % this choice corresponds to ratio cut
            [~, assignment, viols] = cosc(W, vertex_weights, k, ML, CL);
            % Faster, but approximate form: 
            % [cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000);
        catch
            assignment = [];
            disp('Error handled')
        end
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
    end
end

%% CECM
if 1
    name = 'CECM';
    if ~containsMethod(name, titles)
        tic;
        option = struct('init',1,'alpha',1,'rho2',1000,'bal',0,'distance',1);
        noise=0;
        [~,~,BetP,J]=CECM(X,k,matConst,option);
        [~,Idx] = max(BetP,[],2);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = Idx;
        titles{ind+1} = name;
    end
end

%% TVClust variational:
if 1
    name = 'TVClust(variational)';
    if ~containsMethod(name, titles)
        SM = -1 * ones(size(E));
        SM(E == 1) = 1;
        SM(E == -1) = 0;
        
        tic;
        cd '../algorithms/tvclust/'
        [assignment] = TVClust_variational(Checked, SM, X);
        cd '../../experiment/'
        ind = length(dataAll);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
    end
end

%% save everything and possibly plot 
try 
    result = calculateResults(dataAll, assignmentAll, titles, Y);
catch 
end

if verboseOutput 
    save(name1);
    %title = ['Experiment= ' experimentName ' | sample size= ' num2str(size(X,1)) ' | sampling rate= ' num2str(rate) ' | constraint confidence= ' num2str(p) ];
    h = plotExperiments(dataAll, assignmentAll, titles, 4, 3);
    saveas(h, name1, 'fig')
    saveas(h, name1, 'tiff')
    saveas(h, name1, 'png')
    saveas(h, name1, 'pdf')
end 
