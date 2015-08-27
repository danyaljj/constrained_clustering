function [pw, ci, ll, ll_ave] = kmns(dtm, kc, varargin)

% KMNS clusters a group of documents, given by the document-term matrix
% (dtm), using the multinomial model-based k-means algorithm.
% Each document cluster is represented by a multinomial model, whose
% parameters are a vector of word probability distributions.
%
% [pw, ci, ll, ll_ave] = kmns(dtm, kc, [[argID,] value, ...])
%
% Examples:
%	[pw,ci,ll,la] = kmns(dtm,kc)
%	[pw,ci,ll,la] = kmns(dtm,kc,'ml')
%	[pw,ci] = kmns(dtm,kc,'maxi',1,'iprobs',pw0,'map')
%	[pw,ci] = kmns(dtm,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
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
%	'training' *(string) training 'map' or 'ml' (default='map'), which can
%	                     be given without the preceeding argID 'training'
%
% OUTPUT:
%	pw     (matrix) model parameters, size kc by nw, each row a cluster
%	ci     (vector) resulting cluster identities, length nd
%	ll     (matrix) log-likelihood matrix, size nd by kc, each column a
%	                cluster, each row a document
%	ll_ave (vector) a trace of log-likelihood objective values over
%	                iterations note that the objective of map learning
%	                is different from that of ml learning (with addition
%	                of log(prior))
%
% References:
%   A. McCallum and K. Nigam, "A comparison of event models for naive Bayes
%   text classification," AAAI Workshop on Learning for Text Categorization,
%   1998.
%   S. Zhong and J. Ghosh, "A comparative study of generative models for
%   document clustering," SDM Workshop on Clustering High Dimensional Data
%   and Its Applications, 2003
%
% This function can also be used to do classification by doing the following
% two steps (suppose we train on dtm0, test on dtm1):
%	pw0 = kmns(dtm0,kc,'maxi',0,'ipart',classid);
%	[pw1,ci] = kmns(dtm1,kc,'maxi',1,'iprobs',pw0);
% and then compare the resulting ci with the initial partition vector classid
% to get the number of misclassified documents.
%
% See also: MIXMNS, DAMNS, SKMNS, KBERNS, KVMFS


% check input arguments
error(nargchk(2,10,nargin));

% size of document-term matrix
[nd, nw] = size(dtm);

% default input arguments
maxi = 20;
relerr = 1e-4;
training = 'map';
ci = []; pw = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
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

% pre-allocate some memory
ll = zeros(nd,kc);
old_ci = zeros(nd,1);

tmp = ones(nd+1,3);
tmp(1:nd,1) = (1:nd)';
tmp(nd+1,:) = [nd kc 0];
HUGE = 1.0e100;

if isempty(pw)
    % initialize multinomial parameters 
    tmp(1:nd,2) = ci;
	stmp = spconvert(tmp(1:nd,:));
	switch training
	case 'ml',
    	pw = stmp' * dtm;
    	spw = sum(pw,2);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	pw = diag(sparse(spw)) * pw;
	case 'map',
    	pw = stmp' * dtm + 1; % to avoid zero probabilities
    	spw = sum(pw,2);
    	spw = sparse(1 ./ spw);
    	pw = diag(spw) * pw;
	end
end

for k = 1 : maxi
 
	% calculate the log-likelihood matrix
	% ll = dtm * log(pw')
	lpw = - HUGE * ones(kc,nw);
	I = find(pw>0); lpw(I) = log(pw(I));
	ll = dtm * lpw';
    
    % re-partition based on the log-likelihoods (ll)
    [y, ci] = max(ll,[],2);
	switch training
	case 'ml', ll_ave(k) = mean(y);
	case 'map', ll_ave(k) = mean(y) + sum(lpw(:))/nd;
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
        if abs(ll_ave(k)-ll_ave(k-1)) < abs(ll_ave(k-1))*relerr, break; end
    end

    % Estimation of multinomial parameters 
	tmp(1:nd,2) = ci;
    if max(ci) < kc
		stmp = spconvert(tmp);
	else
		stmp = spconvert(tmp(1:nd,:));
	end
	switch training
	case 'ml', % ML estimation
    	pw = stmp' * dtm;
    	spw = sum(pw,2);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	pw = diag(sparse(spw)) * pw;
	case 'map', % MAP estimation
    	pw = stmp' * dtm + 1; % to avoid zero probabilities
    	spw = sum(pw,2);
    	spw = sparse(1 ./ spw);
    	pw = diag(spw) * pw;
	end

end

return;

