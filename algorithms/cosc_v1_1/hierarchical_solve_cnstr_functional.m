function [fold,FctValOuter,totIters] = hierarchical_solve_cnstr_functional(W,Q,fold,deg,gamma,FISTA,maxiterations,verbosity)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

    [fold,FctValOuter,totIters] = solve_cnstr_functional(W,Q,fold,deg,gamma,maxiterations,verbosity);

end

% description of vague variable names:
% rvalold: alpha (in the theorem)


% if(nargin<8)
%     verbosity=3;
% end
% 
% %totIters = 0;
% %assert(isnumeric(W) && issparse(W),'Wrong usage. W should be sparse and numeric.');
% 
% if length(unique(fold)) <= 1
%     fold(ceil(rand(1)*length(fold))) = fold(1)-1; % set one element to zero since we get all ones vector.. this is because the previous solution has all ones on the current
%     % sparsified graph.
% end
% 
% [ix,jx,wval]=find(W);
% 
% W2=triu(W);
% rvalold=zeros(length(ix),1);
% %    rvalold=zeros(length(ix)/2,1);
% pars.MAXITER=100;
% pars.epsilon = 1E-6;
% 
% 
% %maxiterations = 1000;
% 
% counter=0;
% FctValOld=inf;
% FctValOuter=[];
% fnew=fold;
% 
% % We need a better way to form the total variation term for cvx.
% % 	W1=triu(W);
% %     [I,J,V]=find(W1);
% %     V2=[V, - V] ;
% %     D=sparse(length(I),length(fold));
% %     for l=1:length(I)
% %       D(l,[I(l),J(l)])=V2(l,:);
% %     end
% 
% n = size(W,1);
% volQ = sum(sum(Q));
% [qix, qjx, qval] = find(Q);
% 
% u = zeros(n,1);
% v = zeros(n,1);
% %u(1) = 1; v(1)=1;
% [fmax, max_i] = max(fold);
% [fmin, min_i] = min(fold);
% u(max_i) = 1;
% v(min_i) = 1;
% 
% dual_Obj = inf;
% MAX_ITERS_FLAG = true;
% [ix_orig, jx_orig, wval_orig] = find(W);
% degree = sparse(ix_orig, 1, 1, n, 1);
% %degree = ones(n,1);
% 
% %gamma = 0;
% Obj = -1;
% 
% if verbosity >= 2
%     if gamma*volQ == 0
%         fprintf('%8sStarting the unconstrained NCut minimization\n', ' ');
%     else
%         fprintf('%8sStarting the constrained NCut minimization\n', ' ');
%     end
% end
% 
% if verbosity >= 3
%     fprintf('%8sOptimal thresholding of the solution of the inner problem is disabled.\n', ' ');
% end
% n = length(fold); volV = sum(deg);
% nSuccessiveDescentIters = 0;
% while MAX_ITERS_FLAG %(counter<maxiterations)
%     
%     % fnew = fnew/norm(fnew,1);
%     
%     % compute current functional value. This is half the total
%     % variation! because D is constructed from upper triangular W.
%     sval = wval.*abs(fnew(ix)-fnew(jx));% + gamma* volQ * (max(fnew) - min(fnew));
%     
%     % No need to check for normalized flag, as deg is appropriately
%     % initialized.
%     %        Pfnew = fnew - (fnew'*deg/sum(deg));
%     %        FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(deg'*abs(Pfnew));
%     
%     [fsort,sortind]=sort(fnew);
%     sdeg = deg(sortind);
%     if size(sdeg,2)~=1, sdeg = sdeg'; end;
%     rcumvols = flipud(cumsum(flipud(sdeg)));
%     %rcumvols2 = sum(sdeg) - [0; cumsum(sdeg(1:n-1))];
%     %assert( sum(abs(rcumvols-rcumvols2)) == 0 );
%     %assert( max(abs(rcumvols-rcumvols2))/max(rcumvols) <= 1e-8);
%     vec = zeros(n,1);
%     vec(sortind) = 2*sdeg.*(volV - rcumvols - [rcumvols(2:end); 0])/volV;
%     FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(vec'*fnew);  % we cancelled factor 0.5 from numerator and denominator.
%     %assert( abs( fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1) - FctVal ) <= 1e-6);
%     if counter==0
%         display(['        Starting function value (asserted): ', num2str(FctVal)]);
%     end
%     
%     % if functional value has not yet decreased, increase maximum number of
%     % inner iterations
%     %         if(FctVal>=FctValOld || Obj > 0)       
%     
%     if(FctVal>=FctValOld)
%         
%         if ~FISTA
%             % We are solving the inner problem using cvx, so if the
%             % functional does not decrease, then it means the previous
%             % optimum is close to zero; so we break.
%             
%             break;
%             
%         end
%         
%         if(verbosity>=3)
%             disp(['           Functional has not decreased. Old: ',num2str(FctValOld,'%1.16f'),' - New: ',num2str(FctVal,'%1.16f'),'. Increasing number of inner iterations to ', num2str(pars.MAXITER*2)]);
%         end
%         nSuccessiveDescentIters = 0;
%         pars.MAXITER=pars.MAXITER*2;
%         if(pars.MAXITER>=maxiterations)
%             if verbosity >= 3
%                 fprintf('%11sMaximum number of iterations reached... exiting the inner problem\n',' ');
%             end
%             break;
%         end
%         
%         fold=foldback;
%         %FctOld=FctValOld;
%         FctVal=FctValOld;
%         vec = vecOld;
%         
%     else
%         nSuccessiveDescentIters = nSuccessiveDescentIters+1;
%         fold = fnew;
%         FctValOuter=[FctValOuter,FctVal];
%         
%         %             if(abs(FctVal-FctValOld)/FctValOld) < 1e-8
%         %                 break;
%         %             end
%         
%         
%         % for FISTA
%         if FISTA
%             foldback=fold;
%             %FctOld=FctValOld;
%         end
%         
%         FctValOld = FctVal;
%         vecOld = vec; 
%     end
%     
%     if nSuccessiveDescentIters >= 20
%         pars.MAXITER=min(pars.MAXITER*2, maxiterations);
%         if verbosity >= 3
%             display(['           Increasing (doubling) the number of FISTA iterations as we got 20 successive descent steps (because they could be very small): ', num2str(pars.MAXITER)]);
%         end
%     end
%         
%     % Compute the subgradient of the denominator.
%     %        Deg = spdiags(deg,0,size(fold,1),size(fold,1));
%     %        Pfold = fold - (fold'*deg/sum(deg));
%     %        vec = Deg*sign(Pfold) - deg*sum( Deg*sign( Pfold ))/sum(deg);
%     %         vec = sign(fold);
%     %vec = Deg*sign(Pfold) - ones(size(fold,1),1)*( deg'*Deg*sign( Pfold ))/sum(deg);        
%     
%     % Compute the subgradient of R_2(f), without factor 1/2.
%     r2 = 2*gamma*sparse(qix,1,qval.*sign(fold(qix) - fold(qjx)), size(W,1),1);
%     
%     % FISTA call here.
%     if FISTA
%         deg_unit_W = sparse(ix, 1, 1);
%         max_deg = max(deg_unit_W);
%         
%         if gamma*volQ == 0
%             [fnew,rvalold,dual_Obj,niter]=mex_solve_inner_problem(W2,FctVal*full(vec)/2, rvalold,pars.MAXITER,pars.epsilon,4,sqrt(max_deg));
%         else
%             [fnew,rvalold,v,dual_Obj,niter]=mex_solve_cnstr_inner_problem ...
%                 (W2,FctVal*full(vec)/2,rvalold,pars.MAXITER,pars.epsilon,4,full(r2)/2,v,gamma*volQ/2, sqrt(max_deg), degree);
%         end
%         
%         if size(unique(fnew))==1
%             break;
%         end
%         dual_Obj = -dual_Obj;
%         
%         fnew = fnew/norm(fnew);
%         
%         %            Obj = sum(abs(D*fnew))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;
%         Obj = sum(wval.*abs(fnew(ix)-fnew(jx)))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;
%         
%         if isnan(dual_Obj) || isnan(Obj)
%             break;
%             %disp('stop');
%         end
%         
%         g = fnew;
%         
%     else
%         
%         niter = 0;
%         cvx_begin quiet
%         %cvx_precision high;
%         variable g(n,1);
%         %dual variable alph;
%         minimize(sum(abs(D*g))+gamma*volQ*(max(g)-min(g))/2 - g'*r2/2 - FctVal*g'*vec/2);
%         %minimize(sum(abs(D*g)) - 0.5*FctVal*g'*vec);
%         subject to
%         norm(g,2)<=1;
%         cvx_end; cvx_optval;
%         
%         fnew = g;
%         Obj = sum(wval.*abs(g(ix)-g(jx)))+gamma*volQ*(max(g)-min(g))/2 - g'*r2/2 - FctVal*g'*vec/2;
%         dual_Obj = Obj;
%         if abs(dual_Obj) < 1e-12
%             break;
%         end
%         
%     end
%     
%     % [f, c1,c2,c3] = opt_thresh_cnstr_functional(fnew, W, deg', Q,gamma, 1);
%     %         [ff, obj] = opt_thresh_cnstr_functional_old(W, deg, gamma, Q, fnew);
%     %         assert( abs(c1+c2+c3-obj) < 1e-8);
%     %         f = cnstr_opt_thresh(fnew,W,W,true, deg', 1,1:n,1:n,0,0,n, Q, Q, inf, sum(sum(Q>0))/2, 1, 1, true);
%     %         f = fnew>0;
%     %         f = fnew;
%     %         f(fnew>0)=1; f(fnew<=0) =0;
%     %         w = W+0.5*full(sum(sum(Q)))-0.5*Q;
%     %         f= opt_thresh(fnew,w, deg',1);
%     %         if f ~= inf
%     % f=f/norm(f);
%     % fnew = f;
%     %         end
%     if verbosity >= 3 && (counter<100 || rem(counter,10) == 0)
%         disp(['        iter: ', num2str(counter), ' FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, g, 1)), ' FctVal after opt_thres: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' iters: ', num2str(niter)]);
%         %             %disp(['iter: ', num2str(counter), 'FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, g, 1)), 'FctVal after opt_thres: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' duality gap is: ', num2str(Obj - dual_Obj)]);
%         %             %disp(['iter: ', num2str(counter), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' duality gap is: ', num2str(Obj - dual_Obj)]);
%     end
%     
%     totIters = counter;
%     %         if abs(Obj) < 1e-6
%     %             maxiterations = 40000;
%     %             %MAX_ITERS_FLAG = false;
%     %         end
%     
%     %         if Obj > 0
%     %             pars.MAXITER=pars.MAXITER*2;
%     %             disp(['Primal objective is still positive. Increasing number of inner iterations to ', num2str(pars.MAXITER)]);
%     %         end
%     if abs(dual_Obj) < 1e-8
%         MAX_ITERS_FLAG = false;
%     end
%     
%     if counter > maxiterations
%         MAX_ITERS_FLAG = false;
%     end
%     
%     counter=counter+1;
% end
% 
% % compute most recent Functional value
% %sval = wval.*abs(fnew(ix)-fnew(jx));% + gamma* volQ * (max(fnew) - min(fnew));
% %Pfnew = fnew - (fnew'*deg/sum(deg));
% %FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(deg'*abs(Pfnew));
% FctVal = fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1);
% 
% if(FctVal<FctValOld) && length(unique(fnew)) ~= 1
%     fold = fnew;
%     FctValOuter=[FctValOuter,FctVal];
%     if verbosity >= 3
%         disp(['        iter: ', num2str(counter), ' FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew, 1)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' duality gap is: ', num2str(Obj - dual_Obj)]);
%     end
% end
% %    pars
% 
% if verbosity >= 3
%     fprintf('%8sFinal functional value: %f\n', ' ', fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fold, 1));
% end
% end
% 
