%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PART 1: Discrete Distributions and the DP ---
%% Prelim 1.0: bernoulli samples
clear; clf
n = 10;
p = 0.5;
x=0:1:10;
px=binopdf(x, n, p);
ex=binornd(n,p,100,1);
hold on;
subplot(2,1,1);
plot(x, px);
subplot(2,1,2);
hist(ex);

%% Prelim 1.1: generate categorical samples
% The categorical distribution puts a distribution over a discrete number
% of possible outcomes. It's the multidimensional generalization of the
% bernouli distribution
p = [0.2,0.3,0.5];
% R = sum(rand(1) > cumsum(p)) + 1
x = find(mnrnd(1,p))

%% Prelim 1.2: generate multinomial samples
% The multinomial puts a distribution over the number of outcomes for each
% type over the course of n trials.  It's the multidimensional
% generalization of the binomial distribution
% Think: input vector is event probabilities, output vector is event counts
n = 1000;
p = [0.2,0.3,0.5];
X = mnrnd(n,p)

%% Prelim 1.3: generate dirichlet samples
% The dirichlet puts a distribution on the probability vector for a
% multinomial distribution according to a collection of observed
% occurrances of each event.  It's the multidimensional generalization of
% the beta distribution.
% Think: input vector is event counts, output vector is event probabilities
alpha = 1;            % a scale, or concentration parameter
H = [1,1,1];      % any funtion that returns a prob for any event
gamma_samples = gamrnd(alpha*H, 1); % generate n gamma random samples
R = gamma_samples/sum(gamma_samples)% normalize the result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PART 2: The Dirichlet Process ---
%% Prelim 2.1: generate DP marginals of given size
% The dirichlet process puts a distribution over infinite partitions of a
% measurable space E.  We can't instantiate a sample from a DP because it
% would be infinitely long, however we can take advantage of the key
% property that gives it it's name: marginal distributions over the
% unobserved outcomes are dirichlet distributed, with parameters alpha
% (scalar), and H (vector or function).
%
% One run of this process will give you a discrete distribution over
% partitions of some size.  If you ran it enough times and compared all
% size k partitions, you'd find they were dirichlet distributed with
% parameter alpha

clear; clf
alpha = 15;
n = 500;
dp = DirichletProcess(alpha);
dp.DP_run(n);
disp('========================================');
fprintf('num clusters = %d (expected %d)\n', dp.get_k, DirichletProcess.get_k_expected(alpha, n));
disp('========================================');
hist(dp.X, -4:0.1:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PART 3: Gaussian Mixtures ---
%% Prelim 3.1: Mixtures: generate m 1-D samples from each of k gaussians
% Just a demo for visualizing our 1-D mixture model

% Underlying generative model: Gaussian means and variances
G = [-6, 1.5;
    0, 0.5;
    4, 1];
W = [0.25; 0.25; 0.5];
n=200;
X=zeros(n,1);
for i =1:n
    g = G(mnrnd(1,W)==1,:);        % sample parameters according to weights
    X(i) = normrnd(g(1), g(2));    % sample values according to params
end
clf
xrange=min(X):0.1:max(X);
hist(X, 100);

%% Prelim 3.2: Combining gaussian generative mixture with DP prior
% This example shows how the gaussian mixture generated from a DP prior
% can change as a function of alpha and the number of iterations it's
% allowed to run.  
clear all; clf

alpha = 10;      % DP alpha
n_iter = 1000;  % number of draws from the DP
n_samples = 350;% number of samples extracted from the mixture

dpmm=DPMM(alpha);
dpmm.run_open(n_iter, n_samples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- PART 4: Setting up for a sampling: gaussian and DP posteriors
%% Prelim 4.1: sampling gaussian parameters from normal-gamma prior
u0 = 0;     % mean                          (for mean prior)
r0 = 1;     % inverse variance              (for mean prior)
a0 = 10;     % gamma alpha (scale) parameter (for variance prior)
b0 = 0.1;   % gamma beta (shape) parameter  (for variance prior)
p = gamrnd(a0, b0);          % the sampled precision

clf;

% PDF
% c=1;
% xrange = -c:0.1:c;
% m = normpdf(xrange, u0, 1/(p*r0));    % the sampled mean
% plot(xrange, m);

% Histogram
m = normrnd(u0, 1/(p*r0), 1000, 1);    % the sampled mean
hist(m)


%% Prelim 4.2: updating to posterior mean with known variance
clear; clf
% The prior
u_0 = 0;                % prior mean
v_0 = 1;                % prior variance
l_0 = 1/v_0;            % prior precision

% Some hard-coded data
% x = [1 7 2 3 4 2 3];
% n = length(x);

% Some Generated data from another distribution
n = 4;         % number of observations
u_gen = 0.8;    % mean of generated data
v_gen = 0.8;    % variance of generated data
x = normrnd(u_gen * ones(n, 1), v_gen);

% ML estimates for u,v,l | x
u_x = mean(x);
v_x = v_gen/n;%var(x, 1);
l_x = 1/v_x;

% Normal Posterior
uu_post = v_gen/(n*v_0+v_gen)*u_0 + n*v_0/(n*v_0+v_gen)*u_x  % mean of posterior mean
lu_post = l_x + l_0;                              % precision of posterior mean
vu_post = 1/lu_post                                % variance of posterior mean

% Sample from posterior
% u_post =  normrnd(uu_post, vu_post);
% v_post = gamrnd(va_post, vb_post);

% Plot prior & posterior means
clf;
xrange = -3:0.005:3;
hold on;
plot(xrange, normpdf(xrange, u_0, v_0), 'b');           % prior mean
plot(xrange, normpdf(xrange, uu_post, vu_post), 'r');   % posterior mean

%% Prelim 4.3: Gaussian posterior on mean given observed data
% The prior
u_0 = 0;                % prior mean
a_0 = 2;                % prior alpha 
b_0 = 1;                % prior beta
v_0 = gamrnd(a_0, b_0); % prior variance
l_0 = 1/v_0;            % prior precision

% Some hard-coded data
x = [1 7 2 3 4 2 3];
n = length(x);

% Some Generated data from another distribution
% n=3;       % number of observations
% u_gen = 8;  % mean of generated data
% v_gen = 1;% variance of generated data
% x = normrnd(u_0 * ones(n, 1), v_0);

% ML estimates for u,v,l | x
u_x = mean(x);
v_x = var(x, 1);
l_x = 1/v_x;

% Normal Posterior
uu_post = v_0/(v_0+v_x)*u_x + (v_x/(v_0+v_x))*u_0;  % mean of posterior mean
lu_post = n*l_x + l_0;                              % precision of posterior mean
vu_post = 1/lu_post;                                % variance of posterior mean
va_post = a_0 + n/2;                                % alpha of posterior variance
vb_post = b_0 + (n/2)*v_x;                          % beta of posterior variance

% Sample from posterior
% u_post =  normrnd(uu_post, vu_post);
% v_post = gamrnd(va_post, vb_post);

% Plot prior & posterior means
clf;
xrange = -2:0.1:2;
subplot(2,1,1);
hold on;
plot(xrange, normpdf(xrange, u_0, v_0), 'b');           % prior mean
plot(xrange, normpdf(xrange, uu_post, vu_post), 'r');   % posterior mean

% Plot prior & posterior variances
xrange = 0:0.1:80;
subplot(2,1,2);
hold on;
plot(xrange, gampdf(xrange, a_0, b_0), 'b');            % prior variance
plot(xrange, gampdf(xrange, va_post, vb_post), 'r');    % posterior variance


%% Prelim 4.4: Gaussian posterior on precision given observed data


%% Prelim 4.5: Sampling the alpha posterior (in a DP)
clear;
n = 5;
k = 10;

incr = n/500;
A=.1:incr:n;
p=zeros(length(A),1);
for i=1:length(A)
    a = A(i);
    %p(a) = gamma(a)/gamma(n+a) * prod(gamma(N + a/k)/gamma(a/k));
    %p(i) = a^(k-3/2)*exp(-1/(2*a))*gamma(a)/gamma(n+a);
    p(i) = (k-3/2)*log(a) - (1/(2*a)) + gammaln(a) - gammaln(n+a);
end

% generate softmaxed distribution
beta = 1;
bp = exp(1).^(beta*p);
pp=bp/sum(bp);

% plot
clf;
% subplot(2,1,1);
plot(A,pp);
fprintf('ML alpha = %f\n', A(pp==max(pp)));

% draw some samples from it
n_samples = 100;
a_samp = zeros(n_samples,1);
for i=1:n_samples
    ai = sum(rand(1) > cumsum(pp)) + 1;
    a_samp(i) = A(ai);
end
% subplot(2,1,2);
% hist(a_samp);

% np=p/sum(p);
% clf;
% plot(A,np);
% fprintf('ML alpha = %f\n', A(np==min(np)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- PART 5: Comparing DP and EM mixture models ---

%% Demo 1: inferring means and variances of k gaussians using EM
% % 2D Case
% X = zeros(600,2); 
% X(1:200,:) = normrnd(0,1,200,2); 
% X(201:400,:) = normrnd(5,2,200,2); 
% X(401:600,:) = normrnd(-7,3,200,2); 
% [W,M,V,L] = EM_GM(X,3,[],[],1,[]);

% 1D Case
G2 = [-8, 1.0;
       8, 1.0];
W2 = [0.5; 0.5];
G4 = [-20, 1.0;
      -4, 1.0;
       0, 1.0
       9, 1.0];
W4 = [0.2 0.4 0.2 0.2];
n=50;
X=zeros(n,1);
for i =1:n
    %g = G2(mnrnd(1,W2)==1,:);        % sample parameters according to weights
    g = G4(mnrnd(1,W4)==1,:);        % sample parameters according to weights
    X(i) = normrnd(g(1), g(2));    % sample values according to params
end

% Their version
% figure(1);
% [W,M,V,L] = EM_GM(X,4,[],[],1,[])
[W,M,V,L] = EM_GM(X,4,[],[],1,[]);

% My version (more sensitive to kmeans initialization for some reason)
% figure(2);
% [Phi, C, W] = gaussian_EM(X, 2, 10);
% [Phi, C, W] = gaussian_EM(X, 4, 10);

%% Demo 2: inferring means and variances of (some) gaussians using DP-MM
% This is the primary purpose of this demo.  Here we'll generate some data
% from a mixture of gaussians, and then try to recover the means and
% variances by Gibbs sampling a DP-MM


clear; clf
dpmm=DPMM(100);
dpmm.test_MCMC(10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_sampler;
G2 = [-8, 1.0;
       8, 1.0];
W2 = [0.5; 0.5];
G4 = [-10, 1.0;
      -4, 1.0;
       0, 1.0
       9, 1.0];
W4 = [0.2 0.4 0.2 0.2];
n=50;
X=zeros(n,1);
for i =1:n
    %g = G2(mnrnd(1,W2)==1,:);        % sample parameters according to weights
    g = G4(mnrnd(1,W4)==1,:);        % sample parameters according to weights
    X(i) = normrnd(g(1), g(2));    % sample values according to params
end
dpmm=DPMM(10);
post = dpmm.run_MCMC(X, 1000);

alpha_post = cat(1,post.alpha{:});
k_post = cat(1,post.k{:});
figure(2);
hold on;
subplot(2,1,1);
hist(k_post, length(k_post));
title('Posterior on k');
subplot(2,1,2);
hist(alpha_post, length(alpha_post));
title('Posterior on alpha');