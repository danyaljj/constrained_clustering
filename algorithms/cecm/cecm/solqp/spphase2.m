%   spphase2 
%
 lamda=(1.-beta)*lamda;
% if gap <= 5*toler;
%   lamda = lamda/2;
% end;
 go=0;
 gg = Q*x+c;
 XX = spdiags(x,0,n,n);
 AA = A*XX;
 XX = XX*Q*XX;
 %dx = ones(n,1)./x;
%
%  Repeatly solve an ellipsoid constrained QP problem by solving a linear
%  system equation until find a positive solution.
%
 
 while go <= 0,
%   DD = sparse(1:n,1:n,(lamda*dx).*dx,n,n);
%
%   u=[Q+DD A';A sparse(m,m)]\[-(Q*x+c)+(lamda/n)*dx;sparse(m,1)];
%
%   u=[Q+DD A';A sparse(m,m)]\[-gg;sparse(m,1)];
   u=[XX+lamda*speye(n,n) AA';AA sparse(m,m)]\[-x.*gg;sparse(m,1)];
   %u(1:n)=x.*u(1:n);
   xx=x+x.*u(1:n);
   go=min(xx);
   if go > 0,
     ob=xx'*Q*xx/2+c'*xx;
     go = min([go obvalue-ob+eps]);
   end;
   lamda=2*lamda;
   if lamda >= (1+abs(obvalue))/toler,
     %disp('The problem seems unbounded.');
     if ~exist('y')
       y=-u(n+1:n+m);
     end
     return
   end; 
 end;
%
 y=-u(n+1:n+m);
 u=u(1:n);
 nora = min(u);
 if nora < 0,
   nora=-alpha/nora;
 elseif nora == 0,
   nora=alpha;
 else
   nora=inf;
 end
%
 u =  x.*u;
 w1 = u'*Q*u;
 w2 = -u'*gg;
 if w1 > 0,
  nora=min([w2/w1,nora]);
 end;
 if nora == inf,
  ob = -inf;
 else
   x =x+nora*u;
   ob=x'*Q*x/2+c'*x;
 end;
 clear u dx xx DD w1 w2
%
%  This is the Phase 2 procedure called by SPSOLQP. 
