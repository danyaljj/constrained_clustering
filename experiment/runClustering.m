function [] = runClustering(X, Y, k, experimentName, rate, p)
% X 
% Y 
% k 
% experimentName 
% rate
% p 

iter = 1;
name1 = [experimentName '_iter=' num2str(iter)];
name1 = strrep(name1, '.', '_');
if ~exist([name1 '.mat'])
    disp(['skipping = '  experimentName])
    dataAll = {};
    assignmentAll = {};
    titles = {};
    calculationTime = {};
        
    % make E matrix
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
    disp('loading the file .... ')
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

%% cluster : k-means
name = 'K-means';
if ~containsMethod(name, titles)
    disp('K-means ....');
    tic;
    [centroid, pointsInCluster, assignment]= kmeans2(X, k);
    calculationTime{ind+1} = toc;
    ind = length(dataAll);
    dataAll{ind+1} = X;
    assignmentAll{ind+1} = assignment;
    titles{ind+1} = name;
    disp('kmeans done!!')
end

%% dp-means
% disp('RDP-means ....');
% [centroid, pointsInCluster, assignment, clustersSize, objs, pointsAll, centroindsAll] = constrained_v2(X, lambda, E, 2, .001, 'Gaussian');
% ind = length(dataAll);
% dataAll{ind+1} = X;
% assignmentAll{ind+1} = assignment;
% titles{ind+1} = 'RDP-means (with initial assignment to one cluster)';

if 1
    name = 'RDP-means';
    if ~containsMethod(name, titles)
        tic;
        try
            [centroid, pointsInCluster, assignment_out, clustersSize, objs, pointsAll, centroindsAll] = constrained_dpmeans_slow_old(X, lambda, E, 2, .001);
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
            [idx centroids iter LCVQE time] = lcvqe(X, k, C_m);
        catch 
            idx = [];
            disp('Error handled')
        end
        ind = length(dataAll);
        dataAll{ind+1} = X;
        calculationTime{ind+1} = toc;
        assignmentAll{ind+1} = idx;
        titles{ind+1} = name;
        %save('output_lcvqe')
        %pause;
    end
end

%% cluster: dp-means:
if 1
    name = 'DP-means';
    if ~containsMethod(name, titles)
        tic;
        disp('DP-means ....');
        disp(['experiment name  = ' experimentName])
        [centroid, pointsInCluster, assignment, clusterSize]= dpmeans(X, lambda);
        ind = length(dataAll);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
        % save('data')
    end
end

%% constrained k-means
% [centroid, pointsInCluster, assignment]= constrained_kmeans(X, E, k);
% ind = length(dataAll);
% dataAll{ind+1} = X;
% if length(assignment) < size(X,1)
%     maxComp = max(assignment);
%     diffSize = size(X,1) - length(assignment);
%     assignment = [assignment; randint(diffSize, 1, maxComp) + 1];
% end
% assignmentAll{ind+1} = assignment;
% titles{ind+1} = 'Constrained K-means';

% dataAll
% assignmentAll
% titles
% plotExperiments(dataAll, assignmentAll, titles, 2, 3)
% return;

%% MPCKMeans
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


if 1
    name = 'Co1SC-NC';
    if ~containsMethod(name, titles)
        %% Co1SC, 2012: constrained 1-spectral
        disp('Co1SC, 2012: constrained 1-spectral ...');
        %   W:      n x n similarity matrix, where n is the number of data points.
        %           It has to be a sparse symmetric matrix, with zeros on the diagonal.
        %   deg:    n x 1 degree vector (can also be seen as vertex weights)
        %           for Normalized cut it is sum(W,2)
        %           for Ratio cut, it is ones(size(W,1),1)
        %   k:      number of clusters
        %   ML:     m x 2 martix, specifying m must-link pairs
        %   CL:     p x 2 matrix, specifying p cannot-link pairs
        tic;
        try
            W_nonsparse = squareform(pdist(X));
            W = sparse(W_nonsparse);
            vertex_weights = sum(W,2); % this choice corresponds to normalized cut
            [cut, assignment, viols] = cosc(W, vertex_weights, k, ML, CL);
            %[cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000);
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

if 1
    name = 'Co1SC-RC';
    if ~containsMethod(name, titles)
        tic;
        try
            vertex_weights = ones(size(W,1),1); % this choice corresponds to ratio cut
            [cut, assignment, viols] = cosc(W, vertex_weights, k, ML, CL);
            %[cut, clusters, viols] = cosc(W, vertex_weights, k, ML, CL, 2, 2, 0, false, 1000);
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

if 0
    name = 'CECM';
    if ~containsMethod(name, titles)
        %%
        % CECM
        tic;
        option = struct('init',1,'alpha',1,'rho2',1000,'bal',0,'distance',1);
        %nbConst=10;
        noise=0;
        %matConst=eye(n);
        %matConst=addNewConstraints(x,y,matConst,nbConst,noise,0);
        [m,g,BetP,J]=CECM(X,k,matConst,option);
        [C,Idx] = max(BetP,[],2);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = Idx;
        titles{ind+1} = name;
    end
end
% if 0
%     %%
%     % save('tmp')
%     option = struct('init',1,'alpha',1,'rho2',1000,'bal',0,'distance',1);
%     Xi=1;
%     % compute distances between objects
%     D=X*X';
%     N=diag(diag(D))*ones(size(D));
%     DistObj=sqrt(N+N'-2*D);
%     % Generate constraints
%     %E2=addConstraints(y,nbConst);
%     [m,BetP,J,ab]=CEVCLUS(DistObj,k,E2,Xi);
%     [C,Idx] = max(BetP,[],2);
%     dataAll{ind+1} = X;
%     assignmentAll{ind+1} = Idx;
%     titles{ind+1} = 'CEVCLUS';
% end

%% TVClust variational :
if 1
    name = 'TVClust(variational)';
    if ~containsMethod(name, titles)
        SM = -1 * ones(size(E));
        SM(E == 1) = 1;
        SM(E == -1) = 0;
        
        tic;
        cd '../tvclust/'
        [assignment] = TVClust_variational(Checked, SM, X);
        % fileID = fopen('..\tvclust\result2.txt','r');
        % assignment = fscanf(fileID,'%d');
        cd '../toy_experiment/'
        ind = length(dataAll);
        calculationTime{ind+1} = toc;
        dataAll{ind+1} = X;
        assignmentAll{ind+1} = assignment;
        titles{ind+1} = name;
    end
end

% if 0
%     SM = -1 * ones(size(E));
%     SM(E == 1) = 1;
%     %SM(E == -1) = [];
%
%     cd '..\tvclust\'
%     [assignment] = TVClust_variational(Checked, SM, X);
%     % fileID = fopen('..\tvclust\result2.txt','r');
%     % assignment = fscanf(fileID,'%d');
%     cd '..\toy_experiment\'
%     ind = length(dataAll);
%     dataAll{ind+1} = X;
%     assignmentAll{ind+1} = assignment;
%     titles{ind+1} = 'DP+MRF';
% end


%% save everything
try 
    result = calculateResults(dataAll, assignmentAll, titles, Y);
catch 
end

save(name1);

% title = ['Experiment= ' experimentName ' | sample size= ' num2str(size(X,1)) ' | sampling rate= ' num2str(rate) ' | constraint confidence= ' num2str(p) ];
% h = plotExperiments(dataAll, assignmentAll, titles, 4, 5, title);
%
% saveas(h, name1, 'fig')
% saveas(h, name1, 'tiff')
% saveas(h, name1, 'png')
% saveas(h, name1, 'pdf')

