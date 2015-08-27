function [mov, ci, cs, cs_ave] = mixvmfs(dat, kc, varargin)

% MIXVMFS clusters a group of L2-normalized (unit length) vectors,
% using the von Mises-Fisher model-based soft k-means algorithm.
% The parameter kappa is set to be constant across different clusters,
% and increases expoentially over iterations (from 5 to 250 by default).
% Each cluster is represented by a vMF model, whose parameters are
% a unit length mean vector. If you want to use the same kappa across
% iterations, use a midium to high value since low kappa gets stuck very
% quickly in the training process.
%
% [mov, ci, cs, cs_ave] = mixvmfs(dat, kc, [[argID,] value, ...])
%
% Examples:
%	[mov,ci,cs,ca] = mixvmfs(dat,kc)
%	[mov,ci,cs,ca] = mixvmfs(dat,kc,'ml')
%	[mov,ci] = mixvmfs(dat,kc,'maxi',1,'imeans',mov0)
%	[mov,ci] = mixvmfs(dat,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
%
% INPUT:
%	dat     (matrix) input data, size d by n, each column a unit length
%	                 vector (use UNITNORM to make unit length vectors)
%	kc      (scalar) number of clusters
%	Valid argument IDs and corresponding values are:
%	'maxi'      (scalar) maximum number of iterations (default=20)
%	'relerr'    (scalar) relative stopping criterion (default=1e-4)
%	'ipart'     (vector) initial partition vector, length n,
%	                     each entry must be in [1:kc] 
%	'iparas'    (struct) initial model parameters mov0
%	                     mov.mu - initial mean vectors, size d by kc
%	                     mov.mw - initial priors, size 1 by kc
%	'kappa'     (vector) length can be either maxi or 1,
%	                     if 1, the same value used for all iterations
%
% OUTPUT:
%	mov    (struct) model parameters: mov.mu - mean vectors,
%	                size d by kc, each column a cluster; mov.mw - mixture
%	                weights (priors), size 1 by kc
%	ci     (vector) resulting cluster identities, length n
%	cs     (matrix) cosine similarity matrix, size n by kc, each column a
%	                cluster, each row a data point
%	cs_ave (vector) a trace of average cosine similarity values over
%	                iterations
%
% References:
%   S. Zhong and J. Ghosh, "A comparative study of generative models for
%   document clustering," SDM Workshop on Clustering High Dimensional Data
%   and Its Applications, May 2003
%
% See also: DAVMFS, KVMFS, SKVMFS, KMNS, UNITNORM


% check input arguments
error(nargchk(2,8,nargin));

% size of document-term matrix
[d, n] = size(dat);

% default input arguments
maxi = 20;
relerr = 1e-4;
kappa = []; ci = []; mov = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'iparas', i = i+1; mov = varargin{i}; init = init+1;
		case 'kappa', i = i+1; kappa = varargin{i};
		otherwise, disp(sprintf('Ignoring invalid argID %s',varargin{i}));
		end
	else
		disp(sprintf('Ignoring invalid argument #%d',i+2));
	end
	i = i+1;
end

if init > 1
	error('ipart and imeans cannot be used simultaneously');
elseif init == 1
	% check size of ci and mov
	if ~isempty(ci)
		if length(ci) ~= n
        	error(sprintf('The length of ci must be %d!',n));
		end
    	if min(ci) < 1 | max(ci) > kc
        	error(sprintf('The values in ci must be between 1 and %d!',kc));
    	end
    	if size(ci,1) == 1, ci = ci'; end
	end
	if ~isempty(mov)
		if size(mov.mu,2) ~= kc | size(mov.mu,1) ~= d
			error('Incorrect size for mov.mu');
		end
		len = sum(mov.mu.^2,1); 
		if min(len) < 1-relerr | max(len) > 1+relerr
			error('Each column vector in mov.mu must be unit length');
		end
		if size(mov.mw,1) ~= 1 | size(mov.mw,2) ~= kc
			error('Incorrect size for mov.mw');
		end
		if sum(mov.mw) - 1 > eps
			error('mov.mw must sum up to 1');
		end
	end
end

% pre-allocate some memory
cs = zeros(n,kc);
old_ci = zeros(n,1);
postp = zeros(n,kc);

if isempty(mov)
    % initialize vMF parameters 
	if isempty(ci)
		r = 1 + (rand(d,kc) - 0.5) / 10;
		mov.mu = repmat(mean(dat,2),[1,kc]);
		mov.mu = mov.mu .* r; mov.mu = unitnorm(mov.mu,1);
	else
		inds = (1:n)' + n*(ci-1);
		postp(inds) = 1;
		mov.mu = dat * postp;
		mov.mu = unitnorm(mov.mu,1);
	end
	mov.mw = ones(1,kc)/kc;
end

if ~isempty(kappa)
	if length(kappa) == 1
		kappa = kappa * ones(1,maxi);
	elseif length(kappa) < maxi
		error('kappa must be of length 1 or maxi');
	end
else
	if maxi <= 1
		kappa = 100;
	else
		kappa = 400*((1:maxi)/maxi);
	end
end

for k = 1 : maxi

	% calculate the cosine similarity matrix
	cs = full(dat' * mov.mu);
	ll = repmat(log(mov.mw),[n,1]);
	ll = ll + cs * kappa(k);
    
    % re-partition based on the log-likelihoods (cs)
    [y, ci] = max(ll,[],2);
	postp = ll - repmat(y,[1,kc]);
	postp = exp(postp);
	postp = diag(sparse(1./sum(postp,2))) * postp;
	y = sum(cs.*postp,2);
	cs_ave(k) = mean(y);
    
    % stopping criteria
	if k > 1 & min(max(postp,[],2)) > 1-relerr
		if abs(cs_ave(k) - cs_ave(k-1)) < abs(cs_ave(k-1))*relerr
			break;
		end
	end

    % Estimation of vMF parameters 
	mov.mw = sum(postp,1);
	mov.mw = mov.mw / sum(mov.mw);
    mov.mu = dat * postp;
	mov.mu = unitnorm(mov.mu,1);

end

return;

