function [mov,ci,cs,cs_ave,en] = davmfs(dat,kc,varargin)

% DAVMFS implements von Mises-Fisher model-based clustering via deterministic
% annealing. The parameter kappa in a vMF model has the function of an
% inverse temperature. There are two loops in this algorithm: the outer
% loop iterates over increasing kappa and the inner one over number of EM
% iterations for each kappa. Typical kappa values are between 1 and 500.
%
% [mov,ci,cs,cs_ave,en] = davmfs(dat,kc,[[argID,] value,...])
%
% Examples:
%	[mov,ci,cs,ca] = davmfs(dat,kc,'kappa',[1 1.5 2.25 4])
%	[mov,ci] = davmfs(dat,kc,'kappa',1,'maxi',1,'imeans',mov0)
%	[mov,ci] = davmfs(dat,kc,'maxi',10,'relerr',1e-3,'ipart',ci0)
%
% INPUT:
%	dat   (matrix) input data, size d by n, each column a unit length
%	               vector (use UNITNORM to make unit length vectors)
%	kc    (scalar) number of clusters
%	Valid argument IDs and corresponding values are:
%	'kappa'     (vector) inverse temperature schedule, length of kappa
%	                     equals number of iterations in the outer loop
%	                     (default=1.1^(0:60))
%	'maxi'      (scalar) maximum number of iterations for the inner loop
%	                     (default=20)
%	'relerr'    (scalar) relative stopping criterion (default=1e-4)
%	'ipart'     (vector) initial partition vector, length n,
%	                     each entry must be in [1:kc] 
%	'iparas'    (struct) initial model parameters mov0
%	                     mov.mu - initial mean vectors, size d by kc
%	                     mov.mw - initial priors, size 1 by kc
%
% OUTPUT:
%	mov    (struct) model parameters: mov.mu - mean vectors,
%	                size d by kc, each column a cluster
%	                mov.mw - mixture weights (priors), size 1 by kc
%	ci     (vector) resulting cluster identities, length n
%	cs     (matrix) cosine similarity matrix, size n by kc, each column a
%	                cluster, each row a data point
%	cs_ave (vector) a trace of average cosine similarity values
%	                over iterations
%	en     (vector) a trace of average posteriror entropy
%
% References:
%   S. Zhong and J. Ghosh, "A comparative study of generative models for
%   document clustering," SDM Workshop on Clustering High Dimensional Data
%   and Its Applications, May 2003
%
% See also: MIXVMFS, KVMFS, SKVMFS


% check input arguments
error(nargchk(2,8,nargin));

% size of document-term matrix
[d, n] = size(dat);

% default input arguments
kappa = 1.1 .^ (0:60); 
maxi = 10;
relerr = 1e-4;
ci = []; mov = [];

% parse varargin
i = 1; init = 0;
while i<=length(varargin)
	if ischar(varargin{i})
		switch varargin{i}
		case 'kappa', i = i+1; kappa = varargin{i};
		case 'maxi', i = i+1; maxi = varargin{i};
		case 'relerr', i = i+1; relerr = varargin{i};
		case 'ipart', i = i+1; ci = varargin{i}; init = init+1;
		case 'iparas', i = i+1; mov = varargin{i}; init = init+1;
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

for t = 1 : length(kappa)
    for i = 1 : maxi
        % posterior probabilities
        cs = full(dat' * mov.mu);
        dll = repmat(log(mov.mw),[n 1]);
        dll = dll + cs * kappa(t); % deterministic annealing

        [y,ci] = max(dll,[],2);
        dll = dll - repmat(y,[1,kc]);
        postp = exp(dll);
		yp = sum(postp,2);
		postp = diag(sparse(1./yp)) * postp;
        
        % reestimation of priors and mean vectors
        mov.mw = sum(postp);
        mov.mw = mov.mw / sum(mov.mw);
        mov.mu = dat * postp;
        mov.mu = unitnorm(mov.mu,1);
        
        % log-likelihood objective (with \kappa = 1)
        y = y + log(yp);
        ll_ave(i) = mean(y);
        %disp(sprintf('i = %d, log-likelihood = %f',i,ll_ave(i)));
        % stop if the objective value does not change (much)
        if i > 1
            if abs(ll_ave(i)-ll_ave(i-1)) < abs(ll_ave(i-1))*relerr
                break;
            end
        end
    end % i = 1 : maxi
    
    cs_ave(t) = mean(sum(cs.*postp, 2));
    en(t) = entroa(postp);
    if en(t) < relerr
        if abs(cs_ave(t) - cs_ave(t-1)) < abs(cs_ave(t-1)) * relerr
            break;
        end
    end

end % t = 1 : length(kappa)

return;


