%   spphase1 
%
% solve the scaled least squares against two vectors 
%
% [i,j,sa]=find(A);
% AX = sparse(i,j,sa.*x(j),m,n);
% clear i j sa;
% AX'\[ones(n,1) x.*c1];
% 
 dx = ones(n,1)./x(1:n);
 DD = sparse(1:n,1:n,dx.*dx,n,n);
 [DD A';A sparse(m,m)]\[dx sparse(n,1); sparse(m,1) a];
%
 y1=ans(n+1:n+m,1);
 y2=ans(n+1:n+m,2);
 clear dx ans DD;
 w1=(1/ob - a'*y1)/(1/ob^2 - a'*y2);
 w2=1/(1/ob^2 - a'*y2);
 y1=y1-w1*y2;
 y2=-w2*y2;
%
 w1=b'*y1;
 w2=b'*y2;
 y1=y1/(1+w1);
 y2=y2-w2*y1;
 u=[x(1:n).*(-y2'*A)';x(n+1)*(1-y2'*a);w2/(1+w1)];
 v=[x(1:n).*(y1'*A)' ;x(n+1)*(y1'*a)  ; 1/(1+w1)];
%
%  update the dual and the objective lower bound
%
 if min(u-z*v)>=0, 
   y = y2+z*y1; 
   z=b'*y; 
 end;
 clear y1 y2 w1 w2;
%
%  find the descent direction
%  
 u=u-z*v-((ob-z)/(n+2))*ones(n+2,1);
 nora=max(u);
%
%  update the solution along the descent direction
%
 if nora==u(n+1),
   alpha=1.;
 end;
 v=ones(n+2,1)-(alpha/nora)*u;
 x=(x.*v(1:n+1))/v(n+2);
 clear u v
%
%return
%
%  This is the Phase 1 procedure called by SPSOLQP. 
