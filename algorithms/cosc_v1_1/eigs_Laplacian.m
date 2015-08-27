function [v, lambda, ones_C] = eigs_Laplacian( W, b )
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

A = b';
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
A = A*B_inverse_square_root;

opts.disp = 0;
%[vmax, lmax] = eigs(L, 1, 'LA', opts);
lmax = 2*max(sum(W,2))/min(b); % We need an upper bound on the maximum eigenvalue of the generalized Laplcian.
L = lmax*speye(n,n) - L;
%display(['bound on lambda max: ', num2str(lmax)]);  %, ' trace of L: ', num2str(sum(diag(L)))]);

u = rand(size(L,1),1);
converged = false;
lambda_old = 0;

counter = 0;

%call cpp code here...
% [v, lambda] = eig_lcnstrs_mex( L, P, u);

AATranspose = A*A';
totalTime = tic;
while ~converged && counter < 1000
    
    counter = counter+1;
    u = u/norm(u);
    Lu = L*u;
    v = Lu - ((A*Lu)/AATranspose)*A';
    lambda = u'*v;
    
    err = abs((lambda-lambda_old)/lambda_old);
    if err < 1e-20 || abs(lmax - lambda)/lmax < 1e-16
        %if abs(lmax - lambda)/lmax < 1e-16
        converged = true;
        %display(['iter: ', num2str(counter), ' lambda: ', num2str(lambda), ' error: ', num2str(err), ' lambda difference: ', num2str((lmax-lambda)/lambda), ' converged: ', num2str(converged)]);
    end
    
    lambda_old = lambda;
    
    if norm(v) < 1e-8
        lambda = lmax;
        break;
    end
    if rem(counter,100) == 0
        %  display(['iter: ', num2str(counter), ' lambda: ', num2str(lambda), ' error: ', num2str(err), ' lambda difference: ', num2str((lmax-lambda)/lambda), ' avg time per iteration: ', num2str(toc(iterTime)/50)]);
        if toc(totalTime) > 3600
            %display('Time out for eig computation');
            %display(['iter: ', num2str(counter), ' lambda: ', num2str(lambda), ' error: ', num2str(err), ' lambda difference: ', num2str((lmax-lambda)/lambda), ' converged: ', num2str(converged)]);
            break;
        end        
    end
    u = v;
end

% Get back the solution for the original problem!
v = B_inverse_square_root*v;
v = v/norm(v);
lambda = lmax - lambda;
%disp('Eigenproblem end');
%fprintf('converged in %d iterations\n', counter);

start = createClustersGeneral(v,W,true,-1,2,b,true);
if (sum(start)>sum(start==0)) start=1-start; end
start=start/sum(start);
ones_C = start;
end
