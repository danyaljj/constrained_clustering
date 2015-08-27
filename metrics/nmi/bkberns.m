function [pw, ci, ll, ll_ave] = bkberns(dtm, kc, varargin)

% BKBERNS clusters a group of documents, given by the document-term matrix
% (dtm), using the multivariate Bernoulli model-based balanced k-means. Each
% document vector (row) is first converted to a binary vector with each entry
% indicating whether or not a word is in the document. Each document cluster
% is represented by a Bernoulli model, whose parameters are a vector of
% probabilities. The l-th entry of the vector indicating the probability of
% the l-th word being present in the cluster.
%
% [pw, ci, ll, ll_ave] = bkberns(dtm, kc, [[argID,] value, ...])
%
% Examples:
%	[pw,ci,ll,la] = bkberns(dtm,kc)
%	[pw,ci,ll,la] = bkberns(dtm,kc,'ml','minn',100)
%	[pw,ci] = bkberns(dtm,kc,'maxi',1,'iprobs',pw0,'map')
%	[pw,ci] = bkberns(dtm,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
%
% INPUT:
%	dtm     (matrix) document-term matrix, size nd by nw,
%	                 each row a document
%	kc      (scalar) number of clusters
%	Valid argument IDs and corresponding values are:
%	'maxi'      (scalar) maximum number of iterations (default=20)
%	'relerr'    (scalar) relative stopping criterion (default=1e-4)
%	'ipart'     (vector) initial partition vector, length nd,
%	                     each entry must be in [1:kc] 
%	'iprobs'    (matrix) initial model parameters, size kc by nw
%	                     Cannot be used simultaneously with 'ipart'
%	'training'  (string) training 'map' or 'ml' (default='map'), which can
%	                     be given without the preceeding argID 'training'
%	'minn'   	(scalar) minimum number of data points for each cluster,
%	                     (default=0, which means to get completely
%	                     balanced clustering
%
% OUTPUT:
%	pw     (matrix) model parameters, size kc by nw, each row a cluster
%	ci     (vector) resulting cluster identities, length nd
%	ll     (matrix) log-likelihood matrix, size nd by kc, each column a
%	                cluster, each row a document
%	ll_ave (vector) a trace of log-likelihood objective values over
%	                iterations note that the objective of map learning
%	                is different from that of ml learning (with the 
%	                addition of log(prior))
%
% References:
%   S. Zhong and J. Ghosh, "Scalable, balanced model-based clustering," 
%   SIAM Inf. Conf. Data Mining, May 2003
%
% See also: KBERNS, BKMNS, BKVMFS

% check input arguments
error(nargchk(2,12,nargin));

% size of document-term matrix
[nd, nw] = size(dtm);

% default input arguments
maxi = 20;
relerr = 1e-4;
minn = 0;
training = 'map';
ci = []; pw = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'minn', i = i+1; minn = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'iprobs', i = i+1; pw = varargin{i}; init = init+1;
		case 'training', i = i+1; training = varargin{i};
		case {'map','ml'}, training = varargin{i};
		otherwise, disp(sprintf('Ignoring invalid argID %s',varargin{i}));
		end
	else
		disp(sprintf('Ignoring invalid argument #%d',i+2));
	end
	i = i+1;
end

if init > 1
	error('ipart and iprobs cannot be used simultaneously');
elseif init == 1
	% check size of ci and pw
	if ~isempty(ci)
		if length(ci) ~= nd
        	error(sprintf('The length of ci must be %d!',nd));
		end
    	if min(ci) < 1 | max(ci) > kc
        	error(sprintf('The values in ci must be between 1 and %d!',kc));
    	end
    	if size(ci,1) == 1, ci = ci'; end
	end
	if ~isempty(pw)
		if size(pw,1) ~= kc | size(pw,2) ~= nw
			error('Incorrect size for pw');
		end
		if min(pw(:)) < 0 | max(pw(:)) > 1
			error('Values in pw out of range [0,1]');
		end
	end
else
    % random partition
    ci = mod((1:nd)', kc) + 1;
    ci = ci(randperm(nd));
end

% binarize the document-term matrix
mask = sparse(dtm>0);

% pre-allocate some memory
ll = zeros(nd,kc);
old_ci = zeros(nd,1);

tmp = ones(nd+1,3);
tmp(1:nd,1) = (1:nd)';
tmp(nd+1,:) = [nd kc 0];
HUGE = 1.0e100;

if isempty(pw)
    % initialize Bernoulli parameters 
    tmp(1:nd,2) = ci;
	stmp = spconvert(tmp(1:nd,:));
	switch training
	case 'ml',
    	pw = stmp' * mask;
    	spw = sum(stmp,1);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	pw = diag(sparse(spw)) * pw;
	case 'map',
    	pw = stmp' * mask + 1; % to avoid zero probabilities
    	spw = sum(stmp,1) + 2;
    	spw = sparse(1 ./ spw);
    	pw = diag(spw) * pw;
	end
end

for k = 1 : maxi

	% calculate the log-likelihood matrix
	% ll = mask * log(pw') + (1-mask) * log(1-pw')
	lpw = - HUGE * ones(kc,nw); lpwi = lpw;
	I = find(pw>0); lpw(I) = log(pw(I));
	I = find(pw<1); lpwi(I) = log(1-pw(I));
	ll = mask * lpw';
    ll = ll + repmat(sum(lpwi,2)',[nd,1]);
    ll = ll - mask * lpwi';

    % re-partition based on the log-likelihoods (ll)
	% balanced partition
	ci = bpart(ll, minn);
	inds = (1:nd)' + nd*(ci-1);
	y = ll(inds);
    %[y, ci] = max(ll,[],2);
	switch training
	case 'ml'
		ll_ave(k) = mean(y);
	case 'map'
		ll_ave(k) = mean(y) + (sum(lpwi(:))+sum(lpw(:)))/nd;
	end

    % stopping criteria
	% if partition does not change
    if sum(old_ci==ci)==nd
        break;
    else
        old_ci = ci;
    end
    if k > 1
		% if the objective does not change much
        if abs(ll_ave(k)-ll_ave(k-1)) < abs(ll_ave(k-1))*relerr
			break;
		end
    end

    % Estimation of Bernoulli parameters 
	tmp(1:nd,2) = ci;
    if max(ci) < kc
		stmp = spconvert(tmp);
	else
		stmp = spconvert(tmp(1:nd,:));
	end
	switch training
	case 'ml', % ML estimation
    	pw = stmp' * mask;
    	spw = sum(stmp,1);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	pw = diag(sparse(spw)) * pw;
	case 'map', % MAP estimation
    	pw = stmp' * mask + 1; % to avoid zero probabilities
    	spw = sum(stmp,1) + 2;
    	spw = sparse(1 ./ spw);
    	pw = diag(spw) * pw;
	end

end

return;

