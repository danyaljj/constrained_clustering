%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this simple example, we assume there is a mixture of 2 dimensional
% gaussian variables, where mu's and covariances are unknown. But we know that the
% covariances are diagonal and isotropic. Therefore what we don't know 
% are mu and the scalar factor sigma. We use the dirichlet process to model
% the problem and do clustering. We assume it has a conjugate prior for
% mu and sigma. Since the likelihood is gaussian, the conjugate prior 
% should have a normal-gamma distribution. More specifically, sigma has gamma
% distribution and mu has multivariate student-t distribution. The purpose
% of using dirichlet process is that we do not want to specify the number
% of components in the mixture, but instead give a prior over 1 to
% infinite. Then the Gibbs sampler is used to draw sample from the posterior
% distribution based on the observations. See [2] for detail about the
% algorithm.
% 
% most implementation is encapsulated in a DirichMix class (DirichMix.m)
%
% distributable under GPL
% written by Zhiyuan Weng, Nov 26 2011
%
%
% Reference:
% [1] D Fink, "A Compendium of Conjugate Priors"
% [2] R Neal, "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef DirichMixConstrained < handle    
    properties
        data; % N*2
        N;  % number of data
        K;  % number of clusters
        c;  % N*1
        cn; % K*1 
        mu; % K*2
        sig;% K*1
        DIM;
        beta;
        mu0;
        E;
		
		itern = 0; % the current number of iterations

        % hyperparameters
		% p(mu,sig|a,b,beta,mu0)
		a     = 5;           % Gauss-Gamma param
		b     = .9;          % Gauss-Gamma param
        alpha = .1;    % determine the probability of create a new cluster
        p = 0.9;
        q = 0.9;
	end
    methods
        function [mu,sig] = PriorSampleGaussGamma(obj)
			sig  = 1/sqrt(gamrnd(obj.a,1/obj.b));  % the definitions of a,b in matlab are different from here
            mu = DirichMix.MultivarStudentRnd(obj.a,obj.b,obj.beta,obj.mu0);
        end
        function SetDimension(obj, DIM)
            obj.beta  = .05*eye(DIM);   % Gauss-Gamma param  % only affect the variance of mu
            obj.mu0   = zeros(1,DIM);       % Gauss-Gamma param  % only affect the mean of mu 
            obj.DIM = DIM;
        end
        function SetE(obj, E)
            obj.E  = E;
        end
        function InputData(obj,data0)
            obj.data = data0;
            obj.N    = size(data0,1);
            obj.c    = ones(obj.N,1); % one means associate to the first cluster
            obj.K    = 1;  % start from K=1
            obj.cn   = obj.N; % 
            obj.mu   = mean(obj.data);
            obj.sig  = sqrt(trace(cov(obj.data)))/10;
			obj.itern    = 0;
        end        
        function DoIteration(obj,iter)
            % iter: number of iterations to do
            if (~exist('iter','var'))
                iter=1;
            end
            fprintf('[#iter]\t# of clusters\t|\t# of data in each cluster\n');
            while(iter)
                iter = iter - 1;
                % randomize the iteration order
                neworder = randperm(obj.N);
                for n1 = 1:obj.N
                    n = neworder(n1);
                    % update distribution
                    obj.cn(obj.c(n)) = obj.cn(obj.c(n)) - 1;
                    
                    % remove cluster associated with none
                    if obj.cn(obj.c(n))==0
                        % move last one to the removed space
                        oldind = obj.c(n);
                        obj.cn( oldind )  = obj.cn(obj.K);
                        obj.mu( oldind,:) = obj.mu(obj.K,:);
                        obj.sig(oldind )  = obj.sig(obj.K);
                        obj.c(obj.c==obj.K)= oldind;
                        obj.K = obj.K - 1;
                    end
                    
                    % computer categorical distribution param
                    catp = zeros(1,obj.K+1);
                    for k = 1:obj.K
                        %compute friends and enemies in cluster
                        idxs = find(obj.c==k);
                        nfriends = 0;
                        nenemies = 0;
                        for i=1:1:size(idxs,1)
                            if obj.E(idxs(i,1),n) == 1
                                nfriends = nfriends+1;
                            end
                            if obj.E(idxs(i,1),n) == -1
                                nenemies = nenemies+1;
                            end
                        end
                        
                        prod = (obj.p/(1-obj.q))^nfriends*((1-obj.p)/obj.q)^nenemies;   
                        catp(k) = obj.cn(k)/(obj.N-1+obj.alpha)*DirichMix.norm2dpdf(obj.data(n,:),obj.mu(k,:),obj.sig(k))*prod;
                    end
                    
                    % computer partition function
                    mu_    = zeros(1,obj.DIM);           % theoretically, could be any value
                    sig_   = (obj.a-1)/obj.b; % theoretically, could be any value
                    
                    % for the only observation obj.data(n,:)
                    y_     = obj.data(n,:);
                    
                    % posterior hyperparameter
                    a0_    = obj.a+1;
                    mu0_   = ((obj.beta+eye(obj.DIM))\(obj.beta*obj.mu0'+eye(obj.DIM)*y_'))';
                    beta0_ = (obj.beta+eye(obj.DIM))/sig_^2;
                    b0_    = obj.b+.5*(y_-obj.mu0)*obj.beta*(y_-obj.mu0)';
                    
                    % partition function = prior*likelihood/posterior
                    lkhood 	 = DirichMix.norm2dpdf(obj.data(n,:),mu_,sig_);
                    prir 	 = DirichMix.NormalGammaPDF(mu_,sig_,obj.a,obj.b,obj.beta,obj.mu0);
                    postrior = DirichMix.NormalGammaPDF(mu_,sig_,a0_,b0_,beta0_,mu0_);
                    Zp = lkhood*prir/postrior;
                    
                    catp(obj.K+1) = obj.alpha/(obj.N-1+obj.alpha) * Zp; %
                    
                    % sample c for c_n
                    obj.c(n) = DirichMix.CatSample(catp,1);
                    
                    % whether c falls in the new category
                    if obj.c(n)  > obj.K   % Add one
                        obj.K    = obj.K + 1;
                        obj.cn(obj.K) = 1; % the new cluster has one element
                        [ mu_, sig_ ] = PriorSampleGaussGamma(obj); % sample from 'pure' prior
                        obj.mu(obj.K,:) = mu_;
                        obj.sig(obj.K)  = sig_;
                    else
                        obj.cn(obj.c(n)) = obj.cn(obj.c(n)) + 1;
                    end
                end
                
                % update mu and sig for each cluster
                for k = 1:obj.K
                    % for normal-gamma prior, calc the posterior with gaussian likelihood
                    % Namely
                    % mu ~ (multivariate)student-t;
                    % 1/sig^2 ~ gamma;
                    y0 = obj.data(obj.c==k,:);
                    yN = obj.cn(k); % number of observation associate with current cluster
                    
                    % posterior hyperparameter
                    a_    = obj.a+yN*2/2;
                    mu_   = ((obj.beta+yN*eye(obj.DIM))\(obj.beta*obj.mu0'+yN*eye(obj.DIM)*mean(y0,1)'))';
                    b_    = obj.b+.5*trace(...
                        (y0-repmat(mean(y0,1),yN,1))*eye(obj.DIM)*(y0-repmat(mean(y0,1),yN,1))'...
                        ) + .5*(mean(y0,1)-obj.mu0)*obj.beta*(mean(y0,1)-obj.mu0)';
                    
                    % sample sig
                    obj.sig(k)  = 1/sqrt(gamrnd(a_,1/b_)+0.001); % sample from gamma
					
					% another posterior hyperparameter
					beta_ = obj.sig(k)^2\(obj.beta+yN*eye(obj.DIM));
					
					% sample mu
                    obj.mu(k,:)   = DirichMix.MultivarStudentRnd(a_,b_,beta_,mu_); % sample from t
                end
                
                % iteration++
                obj.itern = obj.itern + 1;
                
                fprintf('[%d]\tK=%d\t|\t',obj.itern,obj.K);
                for k = 1:obj.K
                    fprintf('%d\t',obj.cn(k));
                end
                fprintf('\n');
                
            end
        end
        function ClusterKPlot(obj,k)  % for debug
            y0 = obj.data(obj.c==k,:);
            scatter(y0(:,1),y0(:,2),'o');
			hold on;
            scatter(obj.mu(k,1),obj.mu(k,2),'*');
        end
        function PlotData(obj)
            scatter(obj.data(:,1),obj.data(:,2),25,obj.c);
        end
    end
    methods(Static)
        function mu  = MultivarStudentRnd(a,b,beta,mu)
            % degree of freedom = 2*a
            % precision = a/b*beta
            % location = mu
            n = 2*a;
			% y ~ N(0,inv(a/b*beta));  
			sigma = 1/sqrt(a/b*beta(1)); % since beta is diagonal
			y   = randn(1,size(mu,2))*sigma;
			% u ~ chi-square(n)        
			u   = chi2rnd(2*a);
			% see multivariate student distribution on wiki
            mu   = y*sqrt(n/u) + mu;
        end
        function x = NormalGammaPDF(mu,sig,a0,b0,beta0,mu0)
            % a0,b0,beta0,mu0 hyperparameters
			x = 1/(2*pi)*sqrt(det(beta0/sig^2))...
				*exp(-1/2/sig^2*(mu-mu0)*beta0*(mu-mu0)')...
				*gampdf(1/sig^2,a0,1/b0);
        end        
        function lk  = norm2dpdf(x,mu,sig) % symmetric
            lk = 1/2/pi/sig/sig*exp(-.5/sig/sig*sum((x-mu).^2));
        end
        function x   = WhiteNormalSample( mu, sig, nrow )
			% generate nrow*1 random vector with N(mu,sig^2*eye)
			% input
			% mu: L*1   mean
			% sig: 1*1  covariance matrix = sig^2*eye(L)
			% nrow: N
			%
			% output
			% x: N*L random normal vector
			x = repmat(mu',nrow,1) + randn(nrow,numel(mu))*sig;
        end
        function a   = CatSample(p,n)
            % draw n samples from categorical distribution p
            m = numel(p);
            z = rand(n,1);
            p = p/sum(p);
            t = cumsum(p);
            a = sum(repmat(z,1,m) >= repmat(t,n,1),2)+1;
        end
	end
    
end

