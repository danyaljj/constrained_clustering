function [mom, ci, ll, la, en] = damns(dtm, kc, varargin)

% DAMNS clusters a group of documents, given by the document-term matrix
% (dtm), using the multinomial model-based deterministic annealing algorithm.
% There are two loops in this algorithm: the outer loop iterates over
% annealing temperature and the inner one over number of EM iterations for
% each temperature.
%
% [mom, ci, ll, la, en] = damns(dtm, kc, [[argID,] value, ...])
%
% Examples:
%	[mom,ci,ll,la,en] = damns(dtm,kc,'beta',[1 2 4 8 16])
%	[mom,ci] = damns(dtm,kc,'beta',[1 2],'maxi',1,'iparas',mom0,'map')
%	[mom,ci] = damns(dtm,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
%
% INPUT:
%	dtm     (matrix) document-term matrix, size nd by nw,
%	                 each row a document
%	kc      (scalar) number of clusters
%	Valid argument IDs and corresponding values are:
%	'beta'      (vector) inverse temperature parameter, length equals
%	                     number of iterations in the outer loop
%	                     (default=.5*1.3.^(0:23))
%	'maxi'      (scalar) maximum number of iterations (default=10)
%	'relerr'    (scalar) relative stopping criterion (default=1e-4)
%	'ipart'     (vector) initial partition vector, length nd,
%	                     each entry must be in [1:kc] 
%	'iparas'    (struct) initial model parameters, mom.pw size kc by nw,
%	                     mom.pc size 1 by kc priors.
%	                     Cannot be used simultaneously with 'ipart'
%	'training'  (string) training 'map' or 'ml' (default='map'), which may
%	                     be given without the preceeding argID 'training'
%
% OUTPUT:
%	mom    (struct) model parameters: mom.pw, size kc by nw, each row a
%	                cluster; mom.pc, size 1 by kc, priors
%	ci     (vector) resulting cluster identities, length nd
%	ll     (matrix) log-likelihood matrix, size nd by kc, each column a
%	                cluster, each row a document
%	la     (vector) a trace of log-likelihood objective values over
%	                iterations note that the objective of map learning
%	                is different from that of ml learning (with the
%	                addition of log(prior))
%	en     (vector) a trace of average posteriror entropy
%
% References:
%   S. Zhong and J. Ghosh, "A comparative study of generative models for
%   document clustering," SDM Workshop on Clustering High Dimensional Data
%   and Its Applications, 2003
%
% See also: MIXMNS, KMNS, DABERNS, DAVMFS


% check input arguments
error(nargchk(2,11,nargin));

% size of document-term matrix
[nd, nw] = size(dtm);

% default input arguments
beta = .5 * 1.3.^(0:23);
maxi = 10;
relerr = 1e-4;
training = 'map';
ci = []; mom = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'beta', i = i+1; beta = varargin{i};
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'beta', i = i+1; beta = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'iparas', i = i+1; mom = varargin{i}; init = init+1;
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
	if ~isempty(mom)
		if ~isstruct(mom)
			error('mom must be a struct with fields mom.pw and mom.pc');
		end
		if size(mom.pw,1) ~= kc | size(mom.pw,2) ~= nw
			error('Incorrect size for mom.pw');
		end
		if min(mom.pw(:)) < 0 | max(mom.pw(:)) > 1
			error('Values in mom.pw out of range [0,1]');
		end
		if size(mom.pc,1) ~= 1 | size(mom.pc,2) ~= kc
			error('Incorrect size for mom.pc');
		end
		if min(mom.pc(:)) < 0 | max(mom.pc(:)) > 1
			error('Values in mom.pc out of range [0,1]');
		end
	end
else
    % random partition
    ci = mod((1:nd)', kc) + 1;
    ci = ci(randperm(nd));
end

% normalize dtm by document lengths
a = sum(dtm,2);
lend = 1/mean(a);
a = 1 ./ a;
da = diag(sparse(a));
dtm = da * dtm;

% pre-allocate some memory
ll = zeros(nd,kc);
old_ci = zeros(nd,1);
postp = zeros(nd,kc);
HUGE = 1.0e100;

if isempty(mom)
    % initialize multinomial parameters 
	inds = (1:nd)' + nd*(ci-1);
	postp(inds) = 1;
	switch training
	case 'ml', % ML estimation
		mom.pc = sum(postp,1);
		mom.pc = mom.pc / sum(mom.pc);
    	mom.pw = postp' * dtm;
    	spw = sum(mom.pw,2);
		J = find(spw>0);
    	spw(J) = 1 ./ spw(J);
    	mom.pw = diag(sparse(spw)) * mom.pw;
	case 'map', % MAP estimation
		mom.pc = sum(postp,1) + 1;
		mom.pc = mom.pc / sum(mom.pc);
    	mom.pw = postp' * dtm + lend; % to avoid zero probabilities
    	spw = sum(mom.pw,2);
    	spw = sparse(1 ./ spw);
    	mom.pw = diag(spw) * mom.pw;
	end
end

for t = 1 : length(beta)

	for k = 1 : maxi
 
		% calculate the log-likelihood matrix
		lpw = - HUGE * ones(kc,nw);
		I = find(mom.pw>0); lpw(I) = log(mom.pw(I));
		ll = dtm * lpw';
		dll = ll*beta(t) + repmat(log(mom.pc),[nd,1]);
    
    	% re-partition based on the log-likelihoods (ll)
    	[y, ci] = max(dll,[],2);
		dll = dll - repmat(y,[1,kc]);
		postp = exp(dll);
		yp = sum(postp,2);

		y = y + log(yp);
		switch training
		case 'ml',
			ll_ave(k) = mean(y);
		case 'map',
			ll_ave(k) = mean(y) + (sum(lpw(:))*lend+sum(log(mom.pc)))/nd;
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

    	% Estimation of multinomial model parameters 
		postp = diag(sparse(1./yp)) * postp;
		switch training
		case 'ml', % ML estimation
			mom.pc = sum(postp,1);
			mom.pc = mom.pc / sum(mom.pc);
    		mom.pw = postp' * dtm;
    		spw = sum(mom.pw,2);
			J = find(spw>0);
    		spw(J) = 1 ./ spw(J);
    		mom.pw = diag(sparse(spw)) * mom.pw;
		case 'map', % MAP estimation
			mom.pc = sum(postp,1) + 1;
			mom.pc = mom.pc / sum(mom.pc);
    		mom.pw = postp' * dtm + lend; % to avoid zero probabilities
    		spw = sum(mom.pw,2);
    		spw = sparse(1 ./ spw);
    		mom.pw = diag(spw) * mom.pw;
		end

	end
    
    switch training
        case 'ml',
            la(t) = mean(sum(ll.*postp,2));
        case 'map',
            la(t) = mean(sum(ll.*postp,2)) + ...
                (sum(lpw(:))*lend+sum(log(mom.pc)))/nd;
    end
	en(t) = entroa(postp);
	if t > 1 & en(t) < relerr
		if abs(la(t) - la(t-1)) < abs(la(t-1)) * relerr
			break;
		end
	end

end % for t = 1 : length(beta)

return;

