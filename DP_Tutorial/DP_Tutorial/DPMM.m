classdef DPMM < handle
    % DPMM - Dirichlet Process Mixture Model
    % 
    % A class for inferring 1D gaussian mixture components from data using
    % a dirichlet process prior. 
    % Cooked up for a Pfunk lab meeting tutorial during a hot week in
    % August 2011
    % 
    % OVERVIEW:
    % Step 1: create an initial distribution over the thetas (n draws from
    % G_{0} using poly-urn process.  Since we haven't observed any data yet,
    % the expected number of clusters is roughly \alpha
    % log(1+\frac{n}{\alpha}).  Note that each of these draws corresponds first
    % to a draw of c_{i} (which identifies the cluster to which the atom
    % belongs), and then, if it belongs to a new cluster, to a draw of the
    % parameters theta_{i} which describe this cluster.
    
    % steps for resampling one of the theta i's:
    % 1: compute the conditional distribution P(theta_{i} | theta_{-i})
    % (theta_{i} removed).
    %   - if we we sample a known cluster, great, just move on
    %   - if we sample a new cluster, then draw it's parameters from the
    %   POSTERIOR distribution based on the prior G_{0} and the datapoint x_{i}
    %   in question (I.E., don't forget to condition on the datapoint)
    %
    % 2: for all clusters, update the theta^{*}_{j} according to the POSTERIOR
    % distribution based on the prior G_{0} and all the data points currently
    % associated with each cluster.  Technically we should probably do this
    % after we resample each datapoint (for the clusters that changed), but
    % it's probably okay to do it once per sweep.
    %
    % 3: Resample alpha from its gamma posterior given the number of data
    % points and components
    %
    % ...
    % 
    % (Repeat until bored)
    % (plot stuff)
    %
    % Author: Jon Scholz
    % Date: 9/1/2011
    
    %% Members
    properties(Access=public)
        % DP Parameters
        alpha = 1;  % concentration parameter
        
        % Base distribution parameters (normal-gamma prior)
        u0 = 0;     % mean                          (for mean prior)
        r0 = 0.1;    % inverse variance              (for mean prior)
        a0 = 10;     % gamma alpha (scale) parameter (for variance prior)
        b0 = 0.1;   % gamma beta (shape) parameter  (for variance prior)
        
        % DP data memebers
        Phi = [];   % Nx2 Matrix of component parameters (mean and precision)
        N = [];     % Nx1 vector of component counts
        
        % MCMC-related data members
        X = [];     % arbitrarily long observed data vector (implies N counts)
        C = [];     % equally long vector of component assignments (implies N counts)
        n_pts = 0;  % number of data points
    end
    
    %% Public methods
    methods(Access=public)
        function obj = DPMM(varargin)
            for i=1:length(varargin)
                obj.parse_arg(varargin{i});
            end
            
            obj.reset;
        end
        
        function reset(obj)
            obj.Phi = [];
            obj.N = [];
            obj.X = [];
            obj.C = [];
            obj.n_pts = 0;
        end
        
        function run_open(obj, n, m)
            % Run the DP for n steps "open loop" to generate mixture
            % components, and then draw m values from this mixture
            obj.reset;
            for i=1:n
                obj.sample_prior;
            end
            x = ones(m,1);
            for j=1:m
                ci =sum(rand(1) > cumsum(obj.N/sum(obj.N))) + 1; % sample a component
                x(j) = obj.sample_component(ci, 1);
            end
            
            obj.X = x;
            obj.n_pts = length(x);
            obj.plot_overlay; % obj.plot_model;
            fprintf('num components = %d\n', length(obj.N(obj.N~=0)));
        end
        
        function post = run_MCMC(obj, data, n_sweeps)
            % Gibbs sample the posterior mixture parameters for n_sweeps
            % given the provided data (in R1).
            % Note: Currently implemented is a 2-phase gibbs sampling
            % approach as described in Neal 2000 (Algorithm 2).  Key idea
            % here is to separate the resampling of component parameters
            % from the resampling of indicator variables.  This lets us
            % resample all indicators before updating the cluster params,
            % thereby decreasing the mixing time.
            rand('seed',0);
            obj.reset;
            obj.n_pts = length(data);
            obj.X = data;
            obj.C = zeros(obj.n_pts,1);
            obj.Phi = zeros(obj.n_pts, 2);  % Initialize param vec to max possible size
            obj.N = zeros(obj.n_pts,1);     % Initialize count vec to max possible size
            
            % set our hyperparameters based on the mean and variance of the
            % data
            obj.u0 = mean(data); %normrnd(mean(data), var(data));
            obj.r0 = 1/var(data); %gamrnd(obj.a0, 1/var(data));
            
            % sample components for our data from the prior
            for i=1:obj.n_pts
                obj.sample_prior;
            end
            
            % Run the thing
            obj.sync_C_to_N;
            post = obj.gibbs_twostep(n_sweeps);
        end
        
        function X = sample_base(obj)
            % Draws one sample from the base distribution H
            [mu, prec] = DPMM.normgamrnd(obj.u0, obj.r0, obj.a0, obj.b0);
            X = [mu, prec];
        end
        
        function sample_prior(obj)
            % Samples component parameters from DP prior, only
            % according to current component proportions and base
            % distribution (see resample_component for version that depends
            % on a data value as well)
            a = obj.alpha;
            n = sum(obj.N);
            
            if rand < a/(n+a)
                % get an available cluster location
                freeidx = obj.get_free_idx;
                
                % draw new sample from H
                obj.Phi(freeidx,:) = obj.sample_base;
                obj.N(freeidx) = 1;
                %fprintf('new cluster at %f\n', obj.Phi(nextidx));
            else
                % draw sample from empirical distribution
                p = obj.N/n;
                ci = sum(rand(1) > cumsum(p)) + 1;
                obj.N(ci) = obj.N(ci) + 1;
                %fprintf('incrementing cluster %f\n', ci);
            end
        end
        
        %% Gibbs sampling functions
        function x = get_likelihood(obj, i)
            % returns an unnormalized likelihood score for the ith data
            % point evaluated at each component in the model
            x = zeros(length(obj.N), 1);
            x(obj.N~=0) = normpdf(obj.X(i), obj.Phi(obj.N~=0,1), obj.Phi(obj.N~=0,2));
        end
        
        function X = sample_base_conditional(obj, x)
            % Draws one sample from the base distribution H, conditioned on
            % the provided observation x.
            u_new = (obj.u0 + x)/2; % since we don't know variance that generated x, just split the difference
            r_new = obj.r0;
            a_new = obj.a0 + 1/2;
            b_new = obj.b0;
%             if sum(isnan([u_new, r_new, a_new, b_new]))>0
%                 disp([u_new, r_new, a_new, b_new]);
%                 error('isnan in normgam params! (sample_base_conditional)');
%             end
            [m, p] = DPMM.normgamrnd(u_new, r_new, a_new, b_new);
%             if p==0
%                 error('got zero precision!');
%             end
            X = [m, p];
        end
        
        function update_component_params(obj, m)
            % Resample the values for component m, conditional on all x
            % currently associated with cm, and the prior H
            
            % skip empty clusters
            if obj.Phi(m,2) == 0 || sum(isnan(obj.Phi(m,:)))>0
                return
            end
            
            xbar = mean(obj.X(obj.C==m)); % mean of cluster
            nm = obj.N(m);          % number of points in cluster
            sm = obj.Phi(m,2);    % precision of cluster
            u_new = (xbar*nm*sm + obj.u0*obj.r0)/(nm*sm+obj.r0);
            r_new = nm*sm+obj.r0;
            a_new = obj.a0 + 1/2;
            b_new = obj.b0;
            
%             if sum(isnan([u_new, r_new, a_new, b_new]))>0
%                 disp([u_new, r_new, a_new, b_new]);
%                 error('isnan in normgam params! (update_component_params)');
%             end
%             if r_new<0
%                 error('got negative precision! (%f)', r_new);
%             end
            
            [mu, prec] = DPMM.normgamrnd(u_new, r_new, a_new, b_new);
            obj.Phi(m, 1) = mu;
            obj.Phi(m, 2) = prec;
        end
        
        function resample_component(obj, i)
            % Sample the component associated with the ith data point,
            % conditional on xi, the component proportions, and the
            % component parameters.
            % This will either assign to another existing cluster, or draw
            % a new one from the base posterior given the ith data point
            
            a = obj.alpha;
            n = obj.n_pts;
            ci = obj.C(i); % cluster index
            
%             if length(obj.Phi) ~= length(obj.N)
%                 error('warning: param and count vectors out of sync!');
%             end
%             if n ~= sum(obj.N)
%                 error('warning: count vector out of sync!');
%             end
%             if length(find(obj.N)) ~= length(unique(obj.C))
%                 error('N and C don''t even agree on the number of represented clusters!');
%             end
%             if obj.N(ci) <= 0
%                 error('count vector thinks this cluster is empty!');
%             end
%             obj.check_NC_sync; % consistency check for obj.N and obj.C
            
            % Decrement sample at this point
            obj.N(ci) = obj.N(ci) - 1;
%             if obj.N(ci) < 0
%                 error('bug!!');
%             end
            
            if obj.N(ci) == 0
                % flag as a dead cluster (so obj.C alignment doesn't shift)
                obj.Phi(ci,:) = NaN;
            end
            
            % sample a component for this point
            if rand < a/(n+a)
                % get an available cluster location
                freeidx = obj.get_free_idx;
                % draw new sample from H|x,c
                obj.Phi(freeidx,:) = obj.sample_base_conditional(obj.X(i)); % append params for a new component
%                 if sum(isnan(obj.Phi(freeidx,:)))>0
%                     error('isnan in Phi params! (update_component_params)');
%                 end
%                 if obj.Phi(freeidx,2)<0
%                     error('got negative precision! (%f)', obj.Phi(freeidx,2));
%                 end
                obj.N(freeidx) = 1; % set the count for this component to 1
                obj.C(i) = freeidx; % associate this data point
                
            else
                % draw sample from empirical distribution | x, X_xc
                %p = obj.N/(n-1);
                p = obj.N .* obj.get_likelihood(i);
                p = p/sum(p);
                
                newidx = sum(rand(1) > cumsum(p)) + 1; % sample a component
                obj.N(newidx) = obj.N(newidx) + 1;  % increment the counter
                obj.C(i) = newidx;                % record the assignment
            end
        end
        
        function resample_alpha(obj)
            % Resamples the main DP alpha parameter, conditional on the
            % current distribution of clusters.  I used a softmax trick
            % here to control how likely we were to sample the ML alpha
            % (turned out to work better)

            n = obj.n_pts;
            k = length(obj.N(obj.N~=0));
                        
            incr = n/500;
            A = incr:incr:n;
            p = zeros(length(A),1);
            for i=1:length(A)
                a = A(i);
                p(i) = (k-3/2)*log(a) - (1/(2*a)) + gammaln(a) - gammaln(n+a);
            end
            
            % now let's softmax this thing and sample from it
            beta = 5; % 1.5
            bp = exp(1).^(beta*p);
            pp=bp/sum(bp);
            
            %clf;
            %plot(A,pp);
            %fprintf('ML alpha = %f\n', A(pp==max(pp)));
            
            ai = sum(rand(1) > cumsum(pp)) + 1;
            obj.alpha = A(ai);
        end
        
        function post = gibbs_twostep(obj, n_sweeps)
            % Implements Neal (2000) algorithm 2: a two-step gibbs sampling
            % procedure for computing the posterior DPMM given a set of
            % observed data.  Phase one resamples all indicators, and phase
            % two resamples component parameters (separating these leads to
            % faster mixing - see neal 2000)
            % Note: this differs from Rasmussen 2000 in that we don't
            % perform the full bayesian update of the base parameters for
            % the mean (u,v) and the variance/precision (a,b).
            
            burn = 100;     % burn in period
            ef = 25;       % sample extract frequency
            
            post.Phi = {};
            post.N = {};
            post.alpha = {};
            post.k = {};
            s=1;
            
            for i=1:n_sweeps
                if mod(i,20)==0
                    fprintf('Sweep %d (alpha = %f, k = %d)\n', i, obj.alpha, obj.get_k);
                end
                
                if i>=burn && mod(i-burn, ef)==0
                    fprintf('--> Extracting sample %d (alpha = %f, k = %d)\n', s, obj.alpha, obj.get_k);
                    obj.plot_model;
                    drawnow
                    post.Phi{s} = obj.Phi;
                    post.N{s} = obj.N;
                    post.alpha{s} = obj.alpha;
                    post.k{s} = obj.get_k;
                    s=s+1;
                end
                
                for j=1:obj.n_pts
                    obj.resample_component(j);
                end
                
                for m=1:length(obj.N)
                    obj.update_component_params(m);
                end
                
                % Resample alpha
                obj.resample_alpha;
            end
        end
        
        %% Utils
        function x = sample(obj)
            % Returns a sample according to current mixture parameters
            % Note: should verify vectors first by calling check_NC_sync
            idx = sum(rand(1) > cumsum(obj.N/obj.n_pts)) + 1; % sample a component
            x = obj.sample_component(idx, 1);
        end
        
        function x = sample_component(obj, m, n)
            % Return n samples from component m
            x = normrnd(obj.Phi(m,1), 1/obj.Phi(m,2), n, 1);
        end
        
        function k = get_k(obj)
            % @return the number of clusters in the current distribution
            k = length(unique(obj.C));
        end
        
        function freeidx = get_free_idx(obj)
            % Book keeping functon: find the first zero entry in count
            % matrix for a new compponent, expanding both vectors as
            % necessary.
            available = find(obj.N==0);
            try
                freeidx = available(1);
            catch err
                freeidx = length(obj.N) + 1;
                % expand parameter vectors as necessary
                obj.N = [obj.N; zeros(20, 1)];      %or DPMM.get_k_expected(obj.alpha, obj.n_pts)
                obj.Phi = [obj.Phi; zeros(20, 2)];  %or DPMM.get_k_expected(obj.alpha, obj.n_pts)
            end
        end
        
        function sync_N_to_C(obj)
            % Synchronize the count vector N with the component assignment
            % vector C...
            cids=unique(obj.C);
            for i=1:length(cids)
                nk = sum(obj.C==cids(i));
                obj.N(cids(i)) = nk;
            end
        end
        
        function sync_C_to_N(obj)
            % Synchronize the component assignment vector C with the count
            % vector N...
            
            cts = [];
            for i=1:length(obj.N)
                for j=1:obj.N(i)
                    cts=[cts; i];
                end
            end
            obj.C = cts(randperm(length(cts)));
        end
        
        function synced = check_NC_sync(obj)
            synced = true;
            %disp('----------------- running consistency check ----------------');
            uc=unique(obj.C);
            for u=1:length(uc)
                expected = obj.N(uc(u));
                actual = sum(obj.C==uc(u));
                if expected ~= actual
                    synced = false;
                    fprintf('expected = %d, actual = %d\n', expected, actual);
                    return
                end
            end
        end
        
        function plot_model(obj)
            clf;
            hold on;
            subplot(2,1,1);
            hist(obj.X, min(100, obj.n_pts));
            dmin = min(obj.X);
            dmax = max(obj.X);
            xrange = dmin:0.1:dmax;
            
            % Plot bumps separately:
            %             for i=1:length(obj.N)
            %                 if obj.N(i)>0
            %                     plot(xrange, 10 * obj.N(i) * normpdf(xrange, obj.Phi(i,1), obj.Phi(i,2)));
            %                 end
            %             end
            % Plot additive mixture:
            y = zeros(1,length(xrange));
            for i=1:length(obj.N)
                if obj.N(i)>0
                    y = y + obj.N(i) * normpdf(xrange, obj.Phi(i,1), obj.Phi(i,2));
                end
            end
            subplot(2,1,2);
            plot(xrange, y);
        end
        
        function plot_overlay(obj)
            % same plot, but overlayed and scaled to fit
            clf;
            hold on;
            
            hist(obj.X, min(100, obj.n_pts));
            dmin = min(obj.X);
            dmax = max(obj.X);
            xrange = dmin:0.1:dmax;
            y = zeros(1,length(xrange));
            for i=1:length(obj.N)
                if obj.N(i)>0
                    y = y + 0.17 * obj.N(i) * normpdf(xrange, obj.Phi(i,1), obj.Phi(i,2)); % 0.17 is arbitrary constant to scale to fit visually!
                end
            end
            plot(xrange, y, 'r');
        end
        
        function test_MCMC(obj, n_sweeps)
            % Generate some data
            G = [-10, 1;
                10, 1];
            W = [0.5; 0.5];
            n=40;
            X=zeros(n,1);
            for i =1:n
                g = G(mnrnd(1,W)==1,:);        % sample parameters according to weights
                X(i) = normrnd(g(1), g(2));    % sample values according to params
            end
            
            % run sampler
            obj.run_MCMC(X, n_sweeps);
        end
    end
    
    %% Private
    methods(Access=private)
        function parse_arg(obj, arg)
            if isscalar(arg)
                obj.alpha = arg;
            else
                obj.Phi = arg;
            end
        end
    end
    
    %% Static
    methods(Access=public, Static)
        function k = get_k_expected(a, n)
            % @return the expected number of clusters given the values of
            % a and n
            k = floor(a * log(1 + n/a));
        end
        
        function [m, p] = normgamrnd(mu, prec, alpha, beta)
            % Sample a mean and precision for this component using the
            % normal-gamma distribution given the specified mean, inverse
            % variance, shape, and scale params
            p = gamrnd(alpha, beta);          % the sampled precision
            m = normrnd(mu, 1/(p*prec));    % the sampled mean
            
            % debugging:
%             if isnan(m)
%                 error('got isnan for m!');
%             end
%             if p<0
%                 error('obtained negative precision! (%f)\n', p);
%             end
%             if isnan(p)
%                 error('got isnan for p!');
%             end
        end
        
        function [m, p] = normgamrnd_vectorized(mu, prec, alpha, beta, n)
            % Sample a mean and precision for this component using the
            % normal-gamma distribution given the specified mean, inverse
            % variance, shape, and scale params
            p = gamrnd(alpha, beta, n, 1);     % the sampled precision
            m = normrnd(mu, 1./(p*prec));   % the sampled mean
            
            % debugging:
%             if isnan(m)
%                 error('got isnan for m!');
%             end
%             if p<0
%                 fprintf('%f %f\n',m, p);
%             end
%             if isnan(p)
%                 error('got isnan for precision!');
%             end
        end
    end
end