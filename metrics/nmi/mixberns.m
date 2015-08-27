function [mob, ci, ll, ll_ave] = mixberns(dtm, kc, varargin)

% MIXBERNS clusters a group of documents, given by the document-term matrix
% (dtm), using the multivariate Bernoulli model-based EM algorithm.
%
% [mob, ci, ll, ll_ave] = mixberns(dtm, kc, [[argID,] value, ...])
%
% Examples:
%	[mob,ci,ll,la] = mixberns(dtm,kc)
%	[mob,ci,ll,la] = mixberns(dtm,kc,'ml')
%	[mob,ci] = mixberns(dtm,kc,'maxi',1,'iparas',mob0,'ml')
%	[mob,ci] = mixberns(dtm,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
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
%	'iparas'    (struct) initial model parameters, mob.pw size kc by nw,
%	                     mob.pc size 1 by kc priors.
%	                     Cannot be used simultaneously with 'ipart'
%	'training' *(string) training 'map' or 'ml' (default='map'), which can
%	                     be given without the preceeding argID 'training'
%
% OUTPUT:
%	mob    (struct) model parameters:
%	                mob.pw, size kc by nw, each row a cluster
%	                mob.pc, size 1 by kc, priors
%	ci     (vector) resulting cluster identities, length nd
%	ll     (matrix) log-likelihood matrix, size nd by kc, each column a
%	                cluster, each row a document
%	ll_ave (vector) a trace of log-likelihood objective values over
%	                iterations note that the objective of map learning
%	                is different from that of ml learning (with the 
%	                addition of log(prior))
%
% References:
%   S. Zhong and J. Ghosh, "A comparative study of generative models for
%   document clustering," SDM Workshop on Clustering High Dimensional Data
%   and Its Applications, 2003
%
% This function can also be used to do classification by doing the following
% two steps (suppose we train on dtm0, test on dtm1):
%	mob0 = mixberns(dtm0,kc,'maxi',0,'ipart',classid);
%	[mob1,ci] = mixberns(dtm1,kc,'maxi',1,'iparas',mob0);
% and then compare the resulting ci with the initial partition vector classid
% to get the number of misclassified documents.
%
% See also: KBERNS, SKBERNS, DABERNS, MIXMNS, MIXVMFS

% check input arguments
error(nargchk(2,10,nargin));

% size of document-term matrix
[nd, nw] = size(dtm);

% default input arguments
maxi = 20;
relerr = 1e-4;
training = 'map';
ci = []; mob = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'iparas', i = i+1; mob = varargin{i}; init = init+1;
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
	error('ipart and iparas cannot be used simultaneously');
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
	if ~isempty(mob)
		if ~isstruct(mob)
			error('mob must be a struct with fields mob.pw and mob.pc');
		end
		if size(mob.pw,1) ~= kc | size(mob.pw,2) ~= nw
			error('Incorrect size for mob.pw');
		end
		if min(mob.pw(:)) < 0 | max(mob.pw(:)) > 1
			error('Values in mob.pw out of range [0,1]');
		end
		if size(mob.pc,1) ~= 1 | size(mob.pc,2) ~= kc
			error('Incorrect size for mob.pc');
		end
		if min(mob.pc(:)) < 0 | max(mob.pc(:)) > 1
			error('Values in mob.pc out of range [0,1]');
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
postp = zeros(nd,kc);
HUGE = 1.0e100;

if isempty(mob)
    % initialize Bernoulli parameters 
	inds = (1:nd)' + nd*(ci-1);
	postp(inds) = 1;
	switch training
	case 'ml', % ML estimation
    	mob.pw = postp' * mask;
    	spw = sum(postp,1);
		mob.pc = spw / sum(spw);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	mob.pw = diag(sparse(spw)) * mob.pw;
	case 'map', % MAP estimation
    	mob.pw = postp' * mask + 1; % to avoid zero probabilities
    	spw = sum(postp,1) + 2;
		mob.pc = spw - 1; mob.pc = mob.pc / sum(mob.pc);
    	spw = sparse(1 ./ spw);
    	mob.pw = diag(spw) * mob.pw;
	end
end

for k = 1 : maxi

	% calculate the log-likelihood matrix
	% ll = mask * log(pw') + (1-mask) * log(1-pw')        --- 'ml'
	%      + sum(log(pw)) + sum(log(1-pw)) + sum(log(pc)) --- 'map'
	lpw = - HUGE * ones(kc,nw); lpwi = lpw;
	I = find(mob.pw>0); lpw(I) = log(mob.pw(I));
	I = find(mob.pw<1); lpwi(I) = log(1-mob.pw(I));
	ll = mask * lpw';
    ll = ll + repmat(sum(lpwi,2)',[nd,1]);
    ll = ll - mask * lpwi';
	ll = ll + repmat(log(mob.pc),[nd,1]);
 
    % re-partition based on the log-likelihoods (ll)
    [y, ci] = max(ll,[],2);
	dll = ll - repmat(y,[1,kc]);
	postp = exp(dll);
	yp = sum(postp,2);
	y = y + log(yp);
	switch training
	case 'ml',
		ll_ave(k) = mean(y);
	case 'map',
		ll_ave(k) = mean(y) + (sum(lpwi(:))+sum(lpw(:))+...
		            sum(log(mob.pc)))/nd;
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

    % Estimation of Bernoulli parameters 
	postp = diag(sparse(1./yp)) * postp;
	switch training
	case 'ml', % ML estimation
    	mob.pw = postp' * mask;
    	spw = sum(postp,1);
		mob.pc = spw / sum(spw);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	mob.pw = diag(sparse(spw)) * mob.pw;
	case 'map', % MAP estimation
    	mob.pw = postp' * mask + 1; % to avoid zero probabilities
    	spw = sum(postp,1) + 2;
		mob.pc = spw - 1;
		mob.pc = mob.pc / sum(mob.pc);
    	spw = sparse(1 ./ spw);
    	mob.pw = diag(spw) * mob.pw;
	end

end

return;


