%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this simple example, we assume there is a mixture of 2 dimensional
% gaussian variables, where mu and covariance are unknown. But we know the
% covariances are diagonal and isotropic. Therefore what we don't know 
% are mu and the scalar factor sigma. We use dirichlet process to model
% the problem and do clustering. We assume it has a conjugate prior for
% mu and sigma. Since the likelihood is gaussian, the conjugate prior 
% should be normal-gamma distribution. More specifically, sigma has gamma
% distribution and mu has multivariate student-t distribution. The purpose
% of using dirichlet process is that we do not want to specify the number
% of components in the mixture, but instead give a prior over 1 to
% infinite. Then Gibbs sampler is used to draw sample from posterior
% distribution based on the observation. See [2] for detail about the
% algorithm.
% 
% most implementation is encapsulated in DirichMix class (DirichMix.m)
%
% distributable under GPL
% written by Zhiyuan Weng, Nov 26 2011
%
%
% Reference:
% [1] D Fink, "A Compendium of Conjugate Priors"
% [2] R Neal, "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% generate sample continues data
size1 = 200; 
mean1 = [2,-1];
cov1 = [1,0.1; 0.1,1];
mean2 = [8,3];
cov2 = [1 .2; 0.2,1];
X = [mvnrnd(mean1, cov1, size1); mvnrnd(mean2, cov2, size1)];
Y = [ones(size1,1)  ; -1*ones(size1,1)]; 
order = randperm(400); 
X = X(order,:); 
Y = Y(order,:); 

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

dirich = DirichMixConstrained; % construct an object of the class
dirich.SetDimension(size(X,2));
dirich.SetE(E);
dirich.InputData(X);
dirich.DoIteration(1000); % 100 iterations
dirich.PlotData



