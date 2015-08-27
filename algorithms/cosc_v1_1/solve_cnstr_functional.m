function [fold,FctValOuter,totIters] = solve_cnstr_functional(W,Q,fold,deg,gamma,maxiterations,verbosity)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

% description of vague variable names:
% rvalold: alpha (in the theorem)

if(nargin<7)
    verbosity=3;
end

if size(W,1)==2 % just two vertex graph
    fold = [0; 1];
    totIters = 0;
    FctValOuter = fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fold, 1);
    return;
end

if length(unique(fold)) <= 1
    fold(ceil(rand(1)*length(fold))) = fold(1)-1; % set one element to zero since we get all ones vector.. this is because the previous solution has all ones on the current
    % sparsified graph.
end

[ix,jx,wval]=find(W);

W2=triu(W);
rvalold=zeros(length(ix)/2,1);
pars.MAXITER=100;
pars.epsilon = 1E-6;

counter=1;
FctValOld=inf;
FctValOuter=[];
fnew=fold;

n = size(W,1);
volQ = sum(sum(Q));
[qix, qjx, qval] = find(Q);


[fmin, min_i] = min(fold);
v = zeros(n,1);
v(min_i) = 1;

dual_Obj = inf;
MAX_ITERS_FLAG = true;
[ix_orig, jx_orig, wval_orig] = find(W);
degree = sparse(ix_orig, 1, 1, n, 1);

Obj = -1;

if verbosity >= 2
    if gamma*volQ == 0
        fprintf('%6sStarting the unconstrained NCut minimization\n', ' ');
    else
        fprintf('%6sStarting the constrained NCut minimization\n', ' ');
    end
end

if verbosity >= 3
    fprintf('%6sOptimal thresholding of the solution of the inner problem is disabled.\n', ' ');
end
n = length(fold); volV = sum(deg);
nSuccessiveDescentIters = 0;
while MAX_ITERS_FLAG %(counter<maxiterations)
    
    % compute current functional value. This is half the total
    % variation! because D is constructed from upper triangular W.
    sval = wval.*abs(fnew(ix)-fnew(jx)); % This is twice the total variation!
    
    % No need to check for normalized flag, as deg is appropriately
    % initialized.
    %         Pfnew = fnew - (fnew'*deg/sum(deg));
    %         FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(deg'*abs(Pfnew));
    
    [fsort,sortind]=sort(fnew);
    sdeg = deg(sortind);
    if size(sdeg,2)~=1, sdeg = sdeg'; end;
    rcumvols = flipud(cumsum(flipud(sdeg)));
    %rcumvols2 = sum(sdeg) - [0; cumsum(sdeg(1:n-1))];
    %assert( max(abs(rcumvols-rcumvols2)) < 1e-8 );
    vec = zeros(n,1);
    vec(sortind) = 2*sdeg.*(volV - rcumvols - [rcumvols(2:end); 0])/volV;
    %assert(sum(vec) < 1e-8);
    
    FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew))) / (vec'*fnew);  % we cancelled factor 0.5 from numerator and denominator.
    
    % if functional value has not yet decreased, increase maximum number of
    % inner iterations
    
    if(FctVal>=FctValOld)
        
        if(verbosity>=3)
            disp(['         Functional has not decreased. Old: ',num2str(FctValOld,'%1.16f'),' - New: ',num2str(FctVal,'%1.16f'),'. Increasing number of inner iterations to ', num2str(pars.MAXITER*2)]);
        end
        
        nSuccessiveDescentIters = 0;
        pars.MAXITER=pars.MAXITER*2;
        if(pars.MAXITER>=maxiterations)
            if verbosity >= 3
                fprintf('%9sMaximum number of iterations reached... exiting the inner problem\n',' ');
            end
            break;
        end
        
        fold=foldback;
        %FctOld=FctValOld;
        FctVal=FctValOld;
        vec = vecOld;
        
    else
        nSuccessiveDescentIters = nSuccessiveDescentIters+1;
        fold = fnew;
        FctValOuter=[FctValOuter,FctVal];
        
        foldback=fold;
        %FctOld=FctValOld;
        FctValOld = FctVal;
        vecOld = vec;
        
    end
    
    if nSuccessiveDescentIters >= 20
        pars.MAXITER=min(pars.MAXITER*2, maxiterations);
        if verbosity >= 3
            display(['           Increasing (doubling) the number of FISTA iterations as we got 20 successive descent steps (because they could be very small): ', num2str(pars.MAXITER)]);
        end
    end
    
    % Compute the subgradient of the denominator.
    %         Deg = spdiags(deg,0,size(fold,1),size(fold,1));
    %         Pfold = fold - (fold'*deg/sum(deg));
    %         vec = Deg*sign(Pfold) - deg*sum( Deg*sign( Pfold ))/sum(deg);
    
    % Compute the subgradient of R_2(f), without factor 1/2.
    r2 = 2*gamma*sparse(qix,1,qval.*sign(fold(qix) - fold(qjx)), size(W,1),1);
    
    % FISTA call here.
    deg_unit_W = sparse(ix, 1, 1);
    max_deg = max(deg_unit_W);
    
    if gamma*volQ == 0
        [fnew,rvalold,dual_Obj,niter]=mex_solve_inner_problem(W2,FctVal*full(vec)/2, rvalold,pars.MAXITER,pars.epsilon,4,sqrt(max_deg));
    else
        [fnew,rvalold,v,dual_Obj,niter]=mex_solve_cnstr_inner_problem ...
            (W2,FctVal*full(vec)/2,rvalold,pars.MAXITER,pars.epsilon,4,full(r2)/2,v,gamma*volQ/2, sqrt(max_deg), degree);
    end
    totIters = counter;
    
    if size(unique(fnew))==1
        break;
    end
    dual_Obj = -dual_Obj;
    
    fnew = fnew/norm(fnew);
    
    %            Obj = sum(abs(D*fnew))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;
    Obj = sum(wval.*abs(fnew(ix)-fnew(jx)))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;
    
    if isnan(dual_Obj) || isnan(Obj) %|| counter == 100
        break;
        %disp('stop');
    end
    
    g = fnew;
    
    if verbosity >= 3 && (counter<100 || rem(counter,10) == 0)
        disp(['        iter: ', num2str(counter), ' FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, g, 1)), ' FctVal after opt_thres: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' iters: ', num2str(niter)]);
    end
    
    if abs(dual_Obj) < 1e-8
        MAX_ITERS_FLAG = false;
    end
    
    if counter > maxiterations
        MAX_ITERS_FLAG = false;
    end
    
    counter=counter+1;
end

% compute most recent Functional value
% 	sval = wval.*abs(fnew(ix)-fnew(jx));
%     Pfnew = fnew - (fnew'*deg/sum(deg));
% 	FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(deg'*abs(Pfnew));

FctVal = fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1);

if(FctVal<FctValOld) && length(unique(fnew)) ~= 1
    fold = fnew;
    FctValOuter=[FctValOuter,FctVal];
    if verbosity >= 3
        disp(['        iter: ', num2str(counter), ' FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' duality gap is: ', num2str(Obj - dual_Obj)]);
    end
end

if verbosity >= 3
    fprintf('%8sFinal functional value: %f\n', ' ', fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fold, 1));
end
end

