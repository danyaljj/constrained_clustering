function [mass,BetP,J,ab]=CEVCLUS(D,c,link,Xi,version)
% CEVCLUS Constrained relational and evidential clustering
% 
% [mass,BetP,J,ab]=CEVCLUS(D,c,link,Xi,version)
% Input
%   D      : dissimilarity matrix, size (n,n)
%   c      : number of clusters
%   link   : constraints, size(nbConst,3). It contains the index of
%            the objects and a value for the type of constraint 
%            (0 = Cannot-Link, 1 = Must-Link)
%   Xi     : tradeoff between the objectif function Jevclus and the 
%            constraints : Jcecm=(1-Xi)Jecm+Xi Jc
%   version: choice of focal elements
%            1 => all the 2^c focal elements
%            2 => singletons + empty set + Omega (default)
%            3 => only singletons
%
% Output
%   mass: mass assignment to each focal element, size (n,f), f=nb 
%         of focal elements
%   BetP: the pignistic probability deduced from mass
%   J   : Stress : stress function
%   ab  : coefficients in the Stress function
%         ab(1)= slope (=a)
%         ab(2)=constant (=b)
%
%  --------------------------------------------------------------------------
%  Author   : Violaine Antoine
%  mail     : violaine.antoine@univ-bpclermont.fr
%  date     : 31-01-2012
%  version  : 2.1
%  reference: V. Antoine, B. Quost, M.-H. Masson and T. Denoeux.
%             CEVCLUS: Evidential clustering with instance-level constraints
%             for proximity data. Soft Computing (to appear), 2013. 
%  --------------------------------------------------------------------------
[n1,n2]=size(D);
if n1~=n2
    error('Matrix D must be a squared matrix');
else
    n=n1;
end;

if nargin==4,
   version=2; % by default, the focals elements are the singletons+ emptyset + Omega
elseif nargin<4,
   disp('ERROR: evclus requires at least 4 parameters')
end;

W=1-eye(n);
I=find(W);
D1=D(I); % the distances are arranged in a vector

if version==3,
   F = eye(c); % focal elements=singletons
elseif version==2,
   F=[eye(c);zeros(1,c);ones(1,c)];
   % focal elements = singletons + empty set + Omega
elseif version==1, % general case, all the 2^c subsets of Omega
    % not recommended unless c very small
   ii=1:2^c;
   F=zeros(length(ii),c);
   for i=1:c, F(:,i)=bitget(ii'-1,i); end;
end;

f=size(F,1);
xi=zeros(f,f);  % the matrix used to compute the degrees of conflict
for i=1:f,
   for j=1:f,
      xi(i,j)=1-max(min(F([i j],:)));
   end;
end;

card=sum(F,2)'; % the cardinality of the subsets 
i0=find(card==0);
if ~isempty(i0), card(i0)=c;end; % the 'cardinality' of the empty set is
% defined as c, this simplifies the calculation of H

% Initialization of credal partition for evlus
alpha=0.1*randn(n,f);

K=conflict(alpha,xi);
K1=K(I); % the degrees of conflict are arranged in a vector, as the distances
ab=[1; 0];


% stress minimization
param=[reshape(alpha,n*f,1);ab(:)];
options=[1 2000 1e-6 20]; % parameters passed to the optimization routine,
%                           see below


% the Stress function
fprintf('EVCLUS initialisation\n');
[P,J,Jev,Jc,it]=harris('coutmds',param,options,D,[],n,Xi,xi,card,W);
param=P;

% the Stress function
fprintf('CEVCLUS\n');
save('tmp2'); 
[P,J,Jev,Jc,it]=harris('coutmds',param,options,D,link,n,Xi,xi,card,W);

alpha=reshape(P(1:n*f),n,f);
[K,mass]=conflict(alpha,xi);
ab=P(n*f+1:end);

BetP=mtobetpF_(mass,F,1);

% ---------------------------------------------------------------------------------
function [x,y,yev,yc,itValid] = harris(fun,x,options,varargin)
%
% Harris optimization method
%
% fun: function to be optimized (usually an M-file fun.m)
% The function 'FUN' should return a scalar function value [F,G]=FUN(x) and 
% the partial derivatives of the function dfun/dx, at the point x.
%
% options(1) = Display parameter : 1 displays some results.
% options(2) = maximum number of iterations
% options(3) = error
% options(4) = number of iterations between two displays
pas = 0.1 * ones(size(x));
a = 1.2;
b = 0.8;
c = 0.5;
ovf = 1e4 ;
unf = 1e-6;

err=[];
it = 0;
gain=1;


[yp,yev,yc,gp]=feval(fun,x,varargin{:});

xp=x;
x = xp - pas .* gp;
itValid=0;
histY=[];
histPRI=[];
while (gain >= options(3)) & (it <= options(2))
  it=it+1;
  [y,yev,yc,g]=feval(fun,x,varargin{:});
  if options(1) >0,
    if rem(it,options(4))==1,
      fprintf(1,'%5.0f %7.4f %6.4g\n',[it y gain]);
    end;
  end;
  if y > yp,    
    x = xp;
    g = gp; 
    pas = pas * c;
    x = x - pas .* g;
    itValid=itValid+1;
  else
    gain=.9*gain+.1*abs(yp-y);
    xp = x;
    test = (g .* gp) >= 0;
    pas = ((test * a) + ((~test) * b)) .* pas;
    pas = (pas<=ovf) .* pas + (pas>ovf) * ovf;
    pas = (pas>=unf) .* pas + (pas<unf) * unf;
    gp = g;
    x = x - pas .* g;
    yp = y;
  end;
end


% ---------------------------------------------------------------------------------
function [J,Jev,Jc,gradJ]=coutmds(param,D,link,n,Xi,xi,card,W);
% 
% The function called by the optimization routine harris
% Computes the objective function J = Stress + lambda*H
% and its derivatives w.r.t alpha(i,k), a and b
% Input:
% param: alpha(i,k), i=1,..,n, k=1,..,f, a and b
% D: dissimialrities
% n: number of objects
% lambda: regularization parameter
% xi: matrix (f,f) used to compute the degrees of conflict
% card: cardinalities of the focal elements
% W : a matrix (n,n) with 0's on the diagonal, and 1's elsewhere
%
% Output:
% J: objective function 
% gradJ: vector (nf+2,1), gradient of J
f=size(xi,1);  % number of focal elements
alpha=reshape(param(1:n*f),n,f);
[K,mass]=conflict(alpha,xi);
ab=param(n*f+1:end);

% Calculation of objective function
if ~isempty(link)
  ptLink=unique(link(:,1:2));
else
  ptLink=[];
end

I=find(W);
K1=K(I);
D1=D(I);
C=1/sum(D1); % the normalizing constant in the Stress function
if min(D1)==0, cst=0.1*mean(D1); else, cst=0; end;
% this is necessary to avoid division by 0 in the Sammon's Stress function
D=D+100*eye(n);

I=C*sum(((ab(1)*K1+ab(2)-D1).^2)./(D1+cst));

if ~isempty(link)
  CLlink=link(find(link(:,3)==0),1:2);
  MLlink=link(find(link(:,3)==1),1:2);
else
  CLlink=[];MLlink=[];
end
 
emptysetPos=find(sum(xi)==f);
if ~isempty(link)
  mixj=[mass(link(:,1),emptysetPos) mass(link(:,2),emptysetPos)];
  xak=diag(sum(~xi)==2);
  Pml=1-(sum(mixj,2)-mixj(:,1).*mixj(:,2))-diag(mass(link(:,1),:)*xak*mass(link(:,2),:)');

  Pcl=diag(mass(link(:,1),:)*(~xi)*mass(link(:,2),:)');
  coefP=link(:,3);
  coefP(coefP==0)=-1;
  P=coefP.*Pml+1-coefP.*Pcl;
else
  P=0;
end

Jev=I;
Jc=mean(P);
J=I+Xi*mean(P);

% Calculation of the gradient
gradalpha=zeros(n,f);
for k=1:n,
	ek=W(k,:).*(ab(1)*K(k,:)+ab(2) - D(k,:));
	ek=ek./(D(k,:)+cst);
	ek(k)=0;
	ek=ek';
	A=mass(k,:)'*mass(k,:);
	v = mass(k,:).*(mass(k,:)-1);
	A=A.*(1-eye(f))+diag(v);	
	B=xi*mass';
	dsk=B'*A;
	dsk=ab(1)*dsk.*repmat(ek,1,f);
	G=sum(dsk,1);
	PartOne = 4*C*G; 

	if ~isempty(link)
	  [pos1 pos2]=find(link(:,1:2)==k);
	  indL=setdiff(unique(link(pos1,1:2)),k);
	  miemptyset=repmat(A(emptysetPos,:),length(indL),1);
	  mjemptyset=repmat(mass(indL,emptysetPos),1,f);
	  cache=repmat((sum(~xi)==2),f,1);
	  gradMLj=miemptyset.*(mjemptyset-1)-((A.*cache)*mass(indL,:)')';
	  gradCLj=mass(indL,:)*(~xi)*A;
	  
	  % modification of the order of the constraints
	  aux=[link(pos1,1:2) link(pos1,3)-1];
	  ordreLink=reshape(aux',1,prod(size(aux)));
	  ordreLink=ordreLink(ordreLink~=k);
	  ordreLink=sortrows(reshape(ordreLink,2,length(ordreLink)/2)');
          coefP=repmat(ordreLink(:,2)+1,1,f);

	  coefP(coefP==0)=-1;

      save tmp3; 
      gradj=coefP.*gradMLj-coefP.*gradCLj;
	else
	  gradj=0;
	end
	gradalpha(k,:)=PartOne+Xi*(sum(gradj,1))/(max(1,size([CLlink; MLlink],1)*2));

end;
grada=2*C*sum((ab(1)*K1+ab(2)-D1).*K1./(D1+cst));
gradb=2*C*sum((ab(1)*K1+ab(2)-D1)./(D1+cst));
gradJ=[reshape(gradalpha,n*f,1);grada;gradb];



%---------------------------------------------------------------------
function [K,mass]=conflict(alpha,xi)
% Calculation of the degrees of conflict for a credal partition alpha
% Input
% alpha: parameters defiing the credal partitions, size (n,f)
% xi: matrix (f,f) for computing the degrees of conflict
%
% Output
% K: matrix (n,n), degrees of conflict
% mass: the bba's, matrix (n,f)

f=size(alpha,2);
mass=exp(-alpha);
mass=mass./repmat(sum(mass,2),1,f);
K=mass*xi*mass';

%---------------------------------------------------------------------
function [out] = mtobetpF_(m, F, flNm)
% computing BetP on ½ from the m vector 
% out = BetP vector: order a,b,c,... only on singletons

[nbVc,nbEf] = size(m);
[nbEf2,nbAt] = size(F);
indVide=(find(sum(F,2)==0));

if nbEf==nbEf2
  out = zeros(nbVc,nbAt);
  
  if m(indVide) > 1-1e-7
    out = ones(nbVc,nbAt)./nbAt; % case of total conflict
				 % => Yager's normalization 
  else
    cardinals = sum(F,2)';
    cardinals(cardinals==0)=1; % to avoid division by 0
    x = m./repmat(cardinals,nbVc,1);	
    
    for i = 1:nbAt
      ind=find(F(:,i)==1);  % ind = list of indices of supersets of atom i
      out(:,i)=sum(x(:,ind),2);
    end
    
    if flNm & ~isempty(indVide)
      out = out./repmat(1-m(:,indVide),1,nbAt);
    end;
    
  end
else
  'ACCIDENT in mtobetpF_: length of the focal set mismatch' 
end