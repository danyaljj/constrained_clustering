function [lambda, v] = eigs_lcnstrs( W, b, A )
% Solves the optimization problem
% minimize/maximize v'*L*v/v'*diag(b)*v subject to A*v = 0; where L is the Laplacian, 
% L = D - W
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
    
    %disp('Eigenproblem start');

    n = size(W,1);
    B=spdiags(b,0,size(W,1),size(W,1));
    D=spdiags(sum(W,2),0,size(W,1),size(W,1));
    L = D-W;
    L = sparse(L);
    
    % First do the linear transformation to convert ellisoid (denominator) to sphere
    % and obtain the equivalent problem
    % minimize w'* B^-1/2*L*B^-1/2 *w/w'*w subject to A*B^-1/2*w = 0, 
    % where w = B^1/2 *v

    b = diag(B);
    B_inverse_square_root = spdiags(b.^-0.5,0,length(B), length(B));
    L = B_inverse_square_root * L * B_inverse_square_root;
    L1 = L;
    L = (L+L')/2;
    %assert( sum(sum( (L-L1).^2 )) <= 1e-12 );
    A = A*B_inverse_square_root;

    opts.disp = 0;
    [vmax, lmax] = eigs(L, 1, 'LA', opts);
    L = lmax*speye(n,n) - L;

    %P = speye(n) - A'*inv(A*A')*A; % Here A is row vector, so A*A' is a scalar.
    %P = speye(n) - A'*((A*A')\A); 
%    P = A'*((A*A')\A);
    
    u = rand(size(L,1),1);
    converged = false;
    lambda_old = 0;
    
%    M = P*L*P;
    %M = P*L;
%    M = L - A'*((A*A')\A)*L; 
    counter = 0;
    
    %call cpp code here...
   % [v, lambda] = eig_lcnstrs_mex( L, P, u);
    while ~converged && counter < 20000
        
        counter = counter+1;
        u = u/norm(u);
        Lu = L*u;
        %v = Lu - P*Lu; 
        v = Lu - A'*((A*A')\(A*Lu)); 
        %v = M*u;
        %v = (P*L*P)\u;
%        v = L*u;
        
        lambda = u'*v;
        
        err = abs((lambda-lambda_old)/lambda_old);
        if err < 1e-32
            converged = true;
        end
        
        lambda_old = lambda;
        
        %if counter==100 && norm(v) < 1e-10
        if norm(v) < 1e-10
            lambda = lmax;
            break;
        end
%         if rem(counter,100) == 0
%             lambda
%         end
        u = v;
    end
    
    % Get back the solution for the original problem!
    v = B_inverse_square_root*v;
    v = v/norm(v);
    lambda = lmax - lambda;
    %disp('Eigenproblem end');
	%fprintf('converged in %d iterations\n', counter);
end

% 
%     v = rand(n,1);
%     v = v/norm(v);
%     v_old = sparse(n,1);
%     beta = 0;
%     
%     converged = false;
%     k=1;
%     while ~converged
%         
%         w = C*v - beta*v_old;
%         alpha = w'*v;
%         w = w - alpha*v;
%         beta = norm(w);
%         v = w/beta;
%    
%         k = k+1;
%         if k>1000
%             converged = true;
%         end
%         
%     end
%    
%     lambda = v'*C*v;
% end
    
    
    
