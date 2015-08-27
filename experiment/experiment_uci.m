clear all; close all; clc;
pathAll('');
rate = 0.03; 
p = 1;  
for iii = 1:5 
    for k_variance = -3:3
    %for rate = [0.01 0.03 0.05]
        %for p = [0.8 0.9 0.95 1]
            %     p = 1;
%                     if rate == 0.01 && iii == 1
%                         continue;   
%                     end

            [X,Y] = readUCIData('iris');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 3;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_iris_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance) ], rate, p);
            end
            %     p = 1;
            [X,Y] = readUCIData('wine');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 3;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_wine_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)], rate, p);                
            end

            %     p = 1;
            [X,Y] = readUCIData('ecoli');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 8;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_eclli_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)], rate, p);
            end

            %     p = 1;
            [X,Y] = readUCIData('glass');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 7;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_glass_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)], rate, p);
            end

            %     p = 1;
            [X,Y] = readUCIData('isomosphere');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 2;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_isomosphere_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance) ], rate, p);
            end

            %     p = 1;
            [X,Y] = readUCIData('balance');
            %[X,Y] = gaussians(200);
            size1 = size(X,1);
            order = randperm(size1);
            X = X(order,:);
            Y = Y(order,:);
            k = 3;
            kprime = k-k_variance;
            if kprime > 0 
                experiment_1(X,Y,kprime, ['Variable_k_balance_p=' num2str(p) '_rate=' num2str(rate) '_i=' num2str(iii) '_kvariance_' num2str(k_variance)], rate, p);
            end
    end
end