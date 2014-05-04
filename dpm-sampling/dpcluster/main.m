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
clear
clc
dirich = DirichMix; % construct an object of the class
K0 = 5;  % number of clusters
Ms = 30; % number of observations associated with each cluster
xp = []; % observations
% generate observations
for k = 1:K0
    [mu,sig] = dirich.PriorSampleGaussGamma; % sample from prior
    xp = [xp;[randn(Ms,2)*sig+repmat(mu,Ms,1),repmat(k,Ms,1)]];
end
subplot(1,2,1)
scatter(xp(:,1),xp(:,2),25,xp(:,3));
title('true clusters');
axis([-5,5,-5,5]);
axis square;
dirich.InputData(xp(:,1:2));
dirich.DoIteration(100); % 100 iterations
subplot(1,2,2)
dirich.PlotData
title('clustering results');
axis([-5,5,-5,5]);
axis square;

