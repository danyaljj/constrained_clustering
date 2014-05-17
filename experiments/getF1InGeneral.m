function [ f1arr ] = getF1InGeneral( X, Y )

k=length(unique(Y));

f1arr= []

%% cluster : k-means 
[centroid, pointsInCluster, assignment]= kmeans2(X, k); 

[p,r] = getbcubed(Y, assignment);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

%% cluster: dp-means: 
% lambda = 10; 
T = mean(X);
[dist, ind] = sort( sqrt(sum((repmat(T,size(X,1),1)-X).^2,2)), 'descend' );

lambda = dist(k); 

[centroid, pointsInCluster, assignment, clusterSize]= dpmeans(X, lambda); 

[p,r] = getbcubed(Y, assignment);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

%% cluster : dpm-variational
T = 50; % maximum number of clusters
[gamma, phi, m, beta, s, p] = variational_dpm(X, 20, T, 1);
[maxVal, clusters] = max(phi);
centers = []; 

for t = 1:T
    xt = X(clusters == t, :);
    if size(xt) ~= 0
        centers = [centers ; m(t,:)];
    end
end

[p,r] = getbcubed(Y, clusters);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

%% cluster : dpm-gibs sampling  : 
dirich = DirichMix; % construct an object of the class
dirich.SetDimension(size(X,2));
dirich.InputData(X);
dirich.DoIteration(150); % 100 iterations
clusters = unique(dirich.c);

[p,r] = getbcubed(Y, dirich.c);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

%% Make E matrix
E = zeros(size(Y, 1), size(Y, 1)); 
Checked = zeros(size(Y, 1), size(Y, 1)); 
randSize = .01 * size(Y, 1) * size(Y, 1); 
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

%% Constrained Gibbs sampling

dirich = DirichMixConstrained; % construct an object of the class
dirich.SetDimension(size(X,2));
dirich.SetE(E);
dirich.InputData(X);
dirich.DoIteration(50); % 100 iterations

clusters = unique(dirich.c);
[p,r] = getbcubed(Y, dirich.c);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

%% constrained bp-means : slow
%transitivity
%if i is is friends with j then all friends of i should be friends with j
%and all friends of j should be friends with i
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

[centroid, pointsInCluster, assignment, clusterSize, objvals, pointsAll, centroindsAll] = constrained_dpmeans_slow(X, lambda, E, 2, 0.001); 
centroid
pointsInCluster
[p,r] = getbcubed(Y,  assignment);
f1arr(end+1) = 2*p*r / (p+r);
f1arr;

end

