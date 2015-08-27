%  [x,y,obhis]=solqp(Q,A,b,c,(toler,beta))
%
%  See the User Guide at the end of this file.
%
%
%  Start Phase 1: try to find an interior feasible point.
%
%
function [x,y,obhis]=solqp(Q,A,b,c,x,verbose,toler,beta)
 if exist('toler') ~= 1 
   toler=1.e-5; 
 end
 if exist('beta') ~= 1 
   beta=0.8;    
 end
 if exist('alpha') ~= 1 
   alpha=0.95;    
 end

 [m,n] = size(A);
 if isempty(x)
   %disp('Search for a feasible point:')
   a=b-A*ones(n,1);
   x=ones(n+1,1);
   z=0;
   ob=x(n+1); 
   obhis=[ob];
   gap = ob - z;
   while gap >= toler,
     spphase1;
     ob=x(n+1);
     obhis=[obhis ob];
     gap = ob - z;
     if z > 0,
       gap = -1;
       disp('The system has no feasible solution.'),
       return
     end;
   end;
   clear a
 else
   ob=0.5*(x'*Q*x)+c'*x;
 end
% 
% Start Phase 2
%
%
 %disp('Search for an optimal solution:');
 alpha = 0.9;
 x=x(1:n);
 comp=rand(n,1);
 [speye(n) A';A sparse(m,m)]\[comp;sparse(m,1)];
 comp=ans(1:n);
 clear ans;
 nora=min(comp./x);
 if nora < 0,
   nora = -.01/nora;
 else
   nora = max(comp./x);
   if nora == 0,
     disp('The problem has a unique feasible point');
     return
   end;
   nora = .01/nora;
 end;
 x = x + nora*comp;
 obvalue=x'*(Q*x)/2+c'*x;
 obhis=[obvalue];
 lower =-inf;
 zhis=[lower];
 gap=1;
 lamda=max([1 abs(obvalue)/sqrt(sqrt(n))]);
 iter=0;

 while gap >= toler,
   iter=iter+1;

   spphase2;
   if ob == -inf,
     gap = 0;
     disp('The problem is unbounded.');
     return
   else 
     obhis=[obhis ob]; 
     comp=Q*x+c-A'*y; 
     if min(comp)>=0
       zhis(iter+1)=ob-x'*comp;
       lower=zhis(iter+1);
       gap=(ob-lower)/(1+abs(ob));
       obvalue=ob;
     else
       zhis(iter+1)=zhis(iter);     
       lower=zhis(iter+1);
       gap=(obvalue-ob)/(1+abs(ob));
       obvalue=ob;       
     end;
   end;
   if iter>200
     fprintf('gap=%f toler=%f',gap,toler);
     end
 end;
 if verbose>2
   disp('A (local) optimal solution is found.');
 end
 return
%
%  This program solves quadratic program in standard form:
%
%     minimize    0.5*(x'*Q*x)+c'*x
%     subject to  A*x=b, x>=0.
%
%  Input 
%      Q: Sparse symmetric objective matrix.
%      A: Sparse constraint left-hand matrix
%      b: constraint right-hand column vector
%      c: objective column vector
%      toler: relative stopping tolerance: the objective value close to 
%             the local optimal one in the range of tolerance. 
%             Default value: 1.e-5.
%      beta : step size: 0 < beta < 1. Default value: .95.
%
%  Output
%     x: (local) optimal solution
%     y: optimal dual solution (Lagrangien multiplier)
%     obhis : objective value history vs iterations
%
%  Subroutines called : spphase1 and spphase2
%
%  To run the program, just type: spsolqp
%
%    This program is the implementation of the interior ellipsoidal trust
%  region and barrier function algorithm with dual solution updating
%  technique in the standard QP form. Two phases are used: the first uses 
%  SPPHASE1 to find an interior feasible point and the second uses SPPHASE2
%  to find a local optimal solution.
%
%  Technical Reference
%  
%  Y. Ye, "An extension of Karmarkar's algorithm and the trust region method
%         for convex quadratic programming," in Progress in Mathematical
%         Programming (N. Megiddo ed.), Springer-Verlag, NY (1989) 49-63.
%
%  Y. Ye, "On affine-scaling algorithm for nonconvex quadratic programming,"
%         Math. Programming 56 (1992) 285-300.
%
%  Comment: Each iteration we solve a linear KKT system like
%
%  ( Q+mu X^{-2}   A^T )(dx) = c'
%  (     A          0  )(dy) = 0
%
%  where X = diag(x)  which is a positive diagonal matrix.



