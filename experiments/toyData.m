clear all; close all; clc;

%% generate sample continues data
size1 = 100;
mean1 = [2,-1];
cov1 = [1,0.1; 0.1,1];
mean2 = [8,3];
cov2 = [1 .2; 0.2,1];
X = [mvnrnd(mean1, cov1, size1); mvnrnd(mean2, cov2, size1)];
Y = [ones(size1,1)  ; -1*ones(size1,1)];
order = randperm(2*size1);
X = X(order,:);
Y = Y(order,:);
k = 2;

%% cluster : k-means
[centroid, pointsInCluster, assignment]= kmeans2(X, k);
Xtmp = X(Y ==1, :);
plot(Xtmp(:, 1), Xtmp(:, 2), 'xr')
hold on;
Xtmp = X(Y ==-1, :);
plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')
for i = 1:k
    plot(centroid(i,1), centroid(i,2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
end


%% cluster: dp-means:
% lambda = 10;
T = mean(X);
[dist, ind] = sort( sqrt(sum((repmat(T,size(X,1),1)-X).^2,2)), 'descend' );

lambda = dist(k);

% figure;
[centroid, pointsInCluster, assignment, clusterSize]= dpmeans(X, lambda);
figure;
% Xtmp = X(Y ==1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xr')
hold on;
% Xtmp = X(Y ==-1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')
for i = 1:clusterSize
    plot(centroid(i,1), centroid(i,2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
    Xtmp = X(assignment ==i, :);
    plot(Xtmp(:, 1), Xtmp(:, 2), 'x',  'color', rand(1,3))
end


%% cluster : dpm
T = 50; % maximum number of clusters
[gamma, phi, m, beta, s, p] = variational_dpm(X, 20, T, 1);
[maxVal, clusters] = max(phi);
centers = [];
figure;
Xtmp = X(Y ==1, :);
plot(Xtmp(:, 1), Xtmp(:, 2), 'xr')
hold on;
Xtmp = X(Y ==-1, :);
plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')

for t = 1:T
    xt = X(clusters == t, :);
    if size(xt) ~= 0
        disp( ['T = ' num2str(t) ' size(xt,1) = ' num2str(size(xt,1)) ' m(t,:) ' num2str(m(t,:)) ])
        centers = [centers ; m(t,:)];
    end
end

for i = 1:size(centers, 1)
    plot(centers(i,1), centers(i,2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
end



%% cluster : dpm-gibs sampling  :
% Daniel: which algorithm is this?
dirich = DirichMix; % construct an object of the class
dirich.SetDimension(size(X,2));
dirich.InputData(X);
dirich.DoIteration(50); % 100 iterations
%dirich.PlotData
clusters = unique(dirich.c);
for i=1:1:size(clusters,1)
    pts = dirich.data(find(dirich.c == clusters(i,1)),:);
    plot(pts(:, 1), pts(:, 2), 'x', 'color', rand(1,3));
    hold on
end

for i=1:1:size(clusters,1)
    pts = dirich.data(find(dirich.c == clusters(i,1)),:);
    m = mean(pts);
    plot(m(:, 1), m(:, 2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10);
    hold on
end

%% Constrained Gibbs sampling

E = zeros(size(Y, 1), size(Y, 1));
Checked = zeros(size(Y, 1), size(Y, 1));
randSize = 0.01 * size(Y, 1) * size(Y, 1);
iterAll = 1;
while(1)
    i1 = randi(size(Y, 1));
    i2 = randi(size(Y, 1));
    if Checked(i1, i2) == 0
        Checked(i1, i2) = 1;
        if Y(i1) == Y(i2)
            E(i1, i2) = 1;
            E(i2, i1) = 1;
        else
            E(i1, i2) = -1;
            E(i2, i1) = -1;
        end
        iterAll = iterAll + 1;
    end
    if( iterAll > randSize)
        break;
    end
end

%transitivity
%if i is is friends with j then all friends of i should be friends with j
%and all friends of j should be friends with i
% for i=1:1:size(E,1)
%     for j=1:1:size(E,2)
%         if( i == j)
%             continue;
%         end
%         if(E(i,j)==1)
%             %             pause;
%             ind1 = find(E(i, :) == 1);
%             ind2 = find(E(:, j) == 1);
%             
%             for k =ind1
%                 if( k == j )
%                     continue;
%                 end
%                 %                 disp(['j = ' num2str(j)  ' k = '  num2str(k) ])
%                 if(  E(j, k)  == -1 || E( k, j)  == -1 )
%                     error('error bitch! ');
%                 end
%                 E(j, k) = 1;
%                 E(k, j) = 1;
%             end
%             for k2 = ind2'
%                 if( k2 == i )
%                     continue;
%                 end
%                 %                 disp(['i = '  num2str(i)  ' k2 = '  num2str(k2) ])
%                 if(  E(i, k2)  == -1 || E( k2, i)  == -1 )
%                     error('error bitch! ');
%                 end
%                 E(i, k2) = 1;
%                 E(k2, i) = 1;
%             end
%         end
%     end
% end
for i=1:1:size(E,1)
    for j=1:1:size(E,1)
        for k=1:1:size(E,1)
            if( i == k || i == j || k == j ) 
                continue; 
            end
            if( E(i, j) == 1 &&  E(j, k) == 1    ) 
                 E(i, k) = 1;  
                 E(k, i) = 1;   
                 if( E(i, k) == -1 ||  E(k, i) == -1) 
                     disp('Bug bitch! 1 ')
                 end 
            end 
            if( E(j, k) == 1 &&  E(k, i) == 1    ) 
                 E(j, i) = 1;
                 E(i, j) = 1;
                 if( E(j, i) == -1 ||  E(i, j) == -1) 
                     disp('Bug bitch! 2 ')
                 end 
            end 
            if( E(k, i) == 1 &&  E(i, j) == 1    ) 
                 E(k, j) = 1;
                 E(j, k) = 1; 
                 if( E(k, j) == -1 ||  E(j, k) == -1) 
                     disp('Bug bitch! 3 ')
                 end                  
            end             
        end
    end
end
%    

dirich = DirichMixConstrained; % construct an object of the class
dirich.SetDimension(size(X,2));
dirich.SetE(E);
dirich.InputData(X);
dirich.DoIteration(50); % 100 iterations

clusters = unique(dirich.c);
for i=1:1:size(clusters,1)
    pts = dirich.data(find(dirich.c == clusters(i,1)),:);
    plot(pts(:, 1), pts(:, 2), 'x', 'color', rand(1,3));
    hold on
end

for i=1:1:size(clusters,1)
    pts = dirich.data(find(dirich.c == clusters(i,1)),:);
    m = mean(pts);
    plot(m(:, 1), m(:, 2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10);
    hold on
end


%% constrained bp-means : fast
E = zeros(size(Y, 1), size(Y, 1));
Checked = zeros(size(Y, 1), size(Y, 1));
randSize = 0.01 * size(Y, 1) * size(Y, 1);
iterAll = 1;
while(1)
    i1 = randi(size(Y, 1));
    i2 = randi(size(Y, 1));
    if Checked(i1, i2) == 0
        Checked(i1, i2) = 1;
        if Y(i1) == Y(i2)
            E(i1, i2) = 1;
            E(i2, i1) = 1;
        else
            E(i1, i2) = -1;
            E(i2, i1) = -1;
        end
        iterAll = iterAll + 1;
    end
    if( iterAll > randSize)
        break;
    end
end

% E = [ 0  1 1 0 ;  1 0 0 0 ;  1 0 0 0 ; 0 0 0 0 ];

%transitivity
%if i is is friends with j then all friends of i should be friends with j
%and all friends of j should be friends with i
for i=1:1:size(E,1)
    for j=1:1:size(E,2)
        if(E(i,j)==1)
            %             pause;
            ind1 = find(E(i, :) == 1);
            ind2 = find(E(:, j) == 1);
            
            for k =ind1
                if( k == j )
                    continue;
                end
                %                 disp(['j = ' num2str(j)  ' k = '  num2str(k) ])
                if E(j, k) == -1 || E(k, j) == -1
                    error('Something is really wrong bitch! ');
                end
                E(j, k) = 1;
                E(k, j) = 1;
            end
            for k2 = ind2'
                if( k2 == i )
                    continue;
                end
                %                 disp(['i = '  num2str(i)  ' k2 = '  num2str(k2) ])
                
                if E(i, k2) == -1 || E(k2, i) == -1
                    error('Something is really wrong bitch! ');
                end
                E(i, k2) = 1;
                E(k2, i) = 1;
            end
        end
    end
end


lambda = 6.1;
xi = 1;
% figure;
[centroid, pointsInCluster, assignment, clusterSize] = constrained_dpmeans_fast(X, lambda, E, xi);
figure;
% Xtmp = X(Y ==1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xr')
hold on;
% Xtmp = X(Y ==-1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')
for i = 1:clusterSize
    plot(centroid(i,1), centroid(i,2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
    Xtmp = X(assignment ==i, :);
    plot(Xtmp(:, 1), Xtmp(:, 2), 'x',  'color', rand(1,3))
end

%% constrained bp-means : slow
E = zeros(size(Y, 1), size(Y, 1));
Checked = zeros(size(Y, 1), size(Y, 1));
randSize = 0.01 * size(Y, 1) * size(Y, 1);
iterAll = 1;
while(1)
    i1 = randi(size(Y, 1));
    i2 = randi(size(Y, 1));
    if i1 == i2
        continue;
    end
    
    if Checked(i1, i2) == 0
        Checked(i1, i2) = 1;
        if Y(i1) == Y(i2)
            E(i1, i2) = 1;
            E(i2, i1) = 1;
        else
            E(i1, i2) = -1;
            E(i2, i1) = -1;
        end
        iterAll = iterAll + 1;
    end
    if( iterAll > randSize)
        break;
    end
end
%transitivity
%if i is is friends with j then all friends of i should be friends with j
%and all friends of j should be friends with i
% for i=1:1:size(E,1)
%     for j=1:1:size(E,2)
%         if(E(i,j)==1)
%             %             pause;
%             ind1 = find(E(i, :) == 1);
%             ind2 = find(E(:, j) == 1);
%             
%             for k =ind1
%                 if( k == j )
%                     continue;
%                 end
%                 %   disp(['j = ' num2str(j)  ' k = '  num2str(k) ])
%                 E(j, k) = 1;
%                 E(k, j) = 1;
%             end
%             for k2 = ind2'
%                 if( k2 == i )
%                     continue;
%                 end
%                 %  disp(['i = '  num2str(i)  ' k2 = '  num2str(k2) ])
%                 E(i, k2) = 1;
%                 E(k2, i) = 1;
%             end
%         end
%     end
% end


for i=1:1:size(E,1)
    for j=1:1:size(E,1)
        for k=1:1:size(E,1)
            if( i == k || i == j || k == j ) 
                continue; 
            end 
            if( E(i, j) == 1 &&  E(j, k) == 1    ) 
                 E(i, k) = 1;  
                 E(k, i) = 1;   
                 if( E(i, k) == -1 ||  E(k, i) == -1) 
                     disp('Bug bitch! 1 ')
                 end 
            end 
            if( E(j, k) == 1 &&  E(k, i) == 1    ) 
                 E(j, i) = 1;
                 E(i, j) = 1;
                 if( E(j, i) == -1 ||  E(i, j) == -1) 
                     disp('Bug bitch! 2 ')
                 end 
            end 
            if( E(k, i) == 1 &&  E(i, j) == 1    ) 
                 E(k, j) = 1;
                 E(j, k) = 1; 
                 if( E(k, j) == -1 ||  E(j, k) == -1) 
                     disp('Bug bitch! 3 ')
                 end                  
            end             
        end
    end
end
%    


lambda = 6.1;
% xi = 1;
% figure;
[centroid, pointsInCluster, assignment, clusterSize, objvals, pointsAll, centroindsAll] = constrained_dpmeans_slow(X, lambda, E, 2, 0.001);
figure(1);
% Xtmp = X(Y ==1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xr')
hold on;
% Xtmp = X(Y ==-1, :);
% plot(Xtmp(:, 1), Xtmp(:, 2), 'xb')
for i = 1:clusterSize
    plot(centroid(i,1), centroid(i,2),'--rs','LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
    Xtmp = X(assignment ==i, :);
    plot(Xtmp(:, 1), Xtmp(:, 2), 'x',  'color', rand(1,3))
end

figure(2);
plot(1:1:size(objvals,2),objvals);

for j = 47:size(centroindsAll,2)
    h = figure(3)
    hold on;
    for i = 1:clusterSize
        centroid1 = centroindsAll{j};
        assignment1  = pointsAll{j};
        plot(centroid1(i,1), centroid1(i,2),'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10)
        Xtmp = X(assignment1 ==i, :);
        plot(Xtmp(:, 1), Xtmp(:, 2), 'x',  'color', rand(1,3))
    end
    pause;
    % hold off
    delete(h);
end
