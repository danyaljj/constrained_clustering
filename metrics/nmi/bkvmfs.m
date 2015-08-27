function [mu, ci, cs, cs_ave] = bkvmfs(dat, kc, varargin)

% KVMFS clusters a group of L2-normalized (unit length) vectors,
% using the von Mises-Fisher model-based balanced k-means algorithm.
% Each cluster is represented by a vMF model, whose parameters are
% a unit length mean vector.
%
% [mu, ci, cs, cs_ave] = bkvmfs(dat, kc, [[argID,] value, ...])
%
% Examples:
%	[mu,ci,cs,ca] = bkvmfs(dat,kc)
%	[mu,ci] = bkvmfs(dat,kc,'maxi',1,'imeans',mu0)
%	[mu,ci] = bkvmfs(dat,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
%
% INPUT:
%	dat     (matrix) input data, size d by n, each column a unit length
%	                 vector (use UNITNORM to make unit length vectors)
%	kc      (scalar) number of clusters
%	Valid argument IDs and corresponding values are:
%	'maxi'      (scalar) maximum number of iterations (default=20)
%	'relerr'    (scalar) relative stopping criterion (default=1e-4)
%	'minn'      (scalar) minimum number of data points in each cluster
%	                     (default=0, which means completely balanced)
%	'ipart'     (vector) initial partition vector, length n,
%	                     each entry must be in [1:kc] 
%	'imeans'    (matrix) initial model parameters, size d by kc
%	                     Cannot be used simultaneously with 'ipart'
%
% OUTPUT:
%	mu     (matrix) model parameters, size d by kc, each column a cluster
%	ci     (vector) resulting cluster identities, length n
%	cs     (matrix) cosine similarity matrix, size n by kc, each column a
%	                cluster, each row a data point
%	cs_ave (vector) a trace of average cosine similarity values over
%	                iterations
%
% References:
%   A. Banerjee and J. Ghosh, "Frequency sensitive competitive learning for
%   clustering on high-dimensional hyperspheres," IJCNN, May 2002
%   S. Zhong and J. Ghosh, "Scalable, balanced model-based clustering,"
%   SIAM Int. Conf. Data Mining, May 2003
%
% See also: KVMFS, BKMNS, BKBERNS


% check input arguments
error(nargchk(2,8,nargin));

% size of document-term matrix
[d, n] = size(dat);

% default input arguments
maxi = 20;
relerr = 1e-4;
minn = 0;
ci = []; mu = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'minn', i = i+1; minn = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'imeans', i = i+1; mu = varargin{i}; init = init+1;
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
	% check size of ci and mu
	if ~isempty(ci)
		if length(ci) ~= n
        	error(sprintf('The length of ci must be %d!',n));
		end
    	if min(ci) < 1 | max(ci) > kc
        	error(sprintf('The values in ci must be between 1 and %d!',kc));
    	end
    	if size(ci,1) == 1, ci = ci'; end
	end
	if ~isempty(mu)
		if size(mu,2) ~= kc | size(mu,1) ~= d
			error('Incorrect size for mu');
		end
		len = sum(mu.^2,1); 
		if min(len) < 1-relerr | max(len) > 1+relerr
			error('Each column vector in mu must be unit length');
		end
	end
end

% pre-allocate some memory
cs = zeros(n,kc);
old_ci = zeros(n,1);

tmp = ones(n+1,3);
tmp(1:n,1) = (1:n)';
tmp(n+1,:) = [n kc 0];

if isempty(mu)
    % initialize vMF parameters 
	if isempty(ci)
		r = 1 + (rand(d,kc) - 0.5) / 10;
		mu = repmat(mean(dat,2),[1,kc]);
		mu = mu .* r; mu = unitnorm(mu,1);
	else
		tmp(1:n,2) = ci;
		stmp = spconvert(tmp(1:n,:));
		mu = dat * stmp; mu = unitnorm(mu,1);
	end
end

for k = 1 : maxi
 
	% calculate the cosine similarity matrix
	cs = full(dat' * mu);
    
    % re-partition based on the log-likelihoods (cs)
	ci = bpart(cs,minn);
	inds = (1:n)' + n*(ci-1);
	y = cs(inds); % [y, ci] = max(cs,[],2);
	cs_ave(k) = mean(y);
    
    % stopping criteria
	% if partition does not change
    if sum(old_ci==ci)==n
        break;
    else
        old_ci = ci;
    end
    if k > 1
		% if the objective does not change much
        if abs(cs_ave(k)-cs_ave(k-1)) < abs(cs_ave(k-1))*relerr, break; end
    end

    % Estimation of vMF parameters 
	tmp(1:n,2) = ci;
    if max(ci) < kc
		stmp = spconvert(tmp);
	else
		stmp = spconvert(tmp(1:n,:));
	end
    mu = dat * stmp;
	mu = unitnorm(mu,1);
end

return;


