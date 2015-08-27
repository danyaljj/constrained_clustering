function [m,g,BetP,J] = CECM(x,K,matConst,parameters);
% CECM Constrained evidential clustering
%
% [m,g,BetP,J] = CECM(x,K,matConst,parameters)
% Input
%   x       : input matrix nxnbAtt
%   K       : number of desired clusters
%   matConst: constraints matrix nxn : 1 values is for ML constraints and -1 
%             for CL contraints. 0 is when there is no constraint. The matrix 
%             become symetric in the algorithm.
%   Optional:
%   parameters.init    : 0=random initialization of the center, 
%                        1=fcm initialization of the center
%   parameters.alpha   : exponent allowing to control the degree of
%                        penalization for the subsets with high cardinality
%   parameters.rho2    : distance of all objects to the empty set
%   parameters.bal     : tradeoff between the objectif function Jecm and the 
%                        constraints : Jcecm=(1-bal)Jecm+bal Jc
%   parameters.distance: 0=euclidean, 1=mahalanobis
%
% Output
%   m   : the mass function
%   g   : the centroids
%   BetP: the pignistic probability deduced from m
%   J   : the objectif function value
%
%  --------------------------------------------------------------------------
%  Author : Violaine Antoine
%  mail   : violaine.antoine@univ-bpclermont.fr
%  date   : 05-26-2010
%  version: 2.1
%  --------------------------------------------------------------------------
[n,nbAtt]=size(x);
if nargin<2
    error('CECM needs two arguments');
end

if nargin==2
  matConst=eye(n);
end

% ------- Check parameters
if ~isfield(parameters,'init') parameters.init=0;end
if ~isfield(parameters,'alpha') parameters.alpha=1;end
if ~isfield(parameters,'rho2') parameters.rho2=100;end
if ~isfield(parameters,'distance') parameters.distance=0;end
if ~isfield(parameters,'bal') parameters.bal=0.5;end

if (parameters.init~=0 & parameters.init~=1) parameters.init=0;end
if parameters.alpha<1 parameters.alpha=1;end
if parameters.rho2<0 parameters.rho2=100;end
if (parameters.distance~=0 & parameters.distance~=1) parameters.distance=0;end
if (parameters.bal<0 | parameters.bal>1) parameters.bal=0.5;end

init=parameters.init;
alpha=parameters.alpha;
rho2=parameters.rho2;
distance=parameters.distance;
bal=parameters.bal;
beta=2;


%--------------- constraint matrix reformulation --------
matContraintes=sign(matConst+matConst'-eye(n));
aux=matConst.*sign(matConst);
aux=max(aux,aux');
matContraintes=matContraintes.*aux;

%--------------- construction of the focal set matrix  ---------
nbFoc=2^K;

k=1:2^K;
F=zeros(length(k),K);
for i=1:K
  F(:,i)=bitget(k'-1,i); 
end


% -- Setting Aeq and beq matrices
Aeq=kron(eye(n),ones(1,nbFoc));
beq=ones(n,1);

% -- Setting lb and ub matrices
lb=zeros(n*nbFoc,1);
ub=ones(n*nbFoc,1);

% -- Setting A and b matrices
A=[];
b=[];

%---------------------- initializations --------------------------------------
% -- centroids initialization
if init==0
  g = rand(K,nbAtt).*repmat(max(x)-min(x),K,1)+repmat(min(x),K,1);
else
  [g res] = fcm(x,K,[NaN,NaN,NaN,0]);
end

% centers calculus for all the subsets
gplus=[];
for i=2:nbFoc
  fi = F(i,:);
  truc = repmat(fi',1,nbAtt);
  gplus = [gplus;sum(g.*truc)./sum(truc)];
end

% compute the euclidean distance
D=[];
for j=1:nbFoc-1
  aux = (x - ones(n,1)*gplus(j,:));
  B = diag(aux*aux');
  D = [D B];
end

% -- compute masses (ECM without constraints)
c = sum(F(2:end,:),2)';
m = zeros(n,nbFoc-1);
for i=1:n
  for j=1:nbFoc-1
    vect1 = D(i,:);
    vect1 = ((D(i,j)*ones(1,nbFoc-1))./vect1).^(1/(beta-1));
    vect2 =  ((c(j)^(alpha/(beta-1)))*ones(1,nbFoc-1))./(c.^(alpha/(beta-1)));
    vect3 = vect1.*vect2 ;
    
    m(i,j)=1/(sum(vect3)+((c(j)^(alpha/(beta-1)))*D(i,j)/rho2)^(1/(beta-1)));
    if isnan(m(i,j))
      m(i,j)=1;
    end
  end
end
m = [abs(ones(n,1)-sum(m,2)) abs(m)];

x0=reshape(m',n*nbFoc,1); % masses used as initialization
[D,S,Smeans]=setDistances(x,F,g,m(:,2:nbFoc),alpha,distance);

% -- Setting f matrix
aux=matContraintes-eye(n);
contraintesML = max(aux,zeros(n));
nbContParObjet = sum(contraintesML,2);
fvide=kron(nbContParObjet,[1; zeros(nbFoc-1,1)]);
%f=ones(n*nbFoc,1).*x0+fvide;
f=fvide;


% -- Setting constraints matrix
nbML=length(find(tril(matContraintes,-1)==1));
nbCL=length(find(tril(matContraintes,-1)==-1));
  
if nbML==0 nbML=1;end
if nbCL==0 nbCL=1;end

MLMat=(sign((F*ones(K,1)-1).^2)-1).^2; % getting focal elements
MLMat=(MLMat*MLMat').*(F*F');
CLMat=sign(F*F');
MLMat=MLMat*-sign(bal)/(2*nbML);
CLMat=CLMat*sign(bal)/(2*nbCL);

% Constraints matrix with the respect of the constraints given in parameters
aux=tril(matContraintes,-1);
contraintesML = max(aux,zeros(n)); % split of the ML and CL constraints
contraintesCL = abs(min(aux,zeros(n)));

contraintesML = sparse(contraintesML);
contraintesCL = sparse(contraintesCL);
  
MLaux=kron(contraintesML,ones(nbFoc));
CLaux=kron(contraintesCL,ones(nbFoc));

contraintesMat=repmat(MLMat,n).*MLaux+repmat(CLMat,n).*CLaux;
contraintesMat=contraintesMat+contraintesMat';


% -- Setting H matrix
aux=D*[zeros(nbFoc-1,1) eye(nbFoc-1)]+[ones(n,1)*rho2 zeros(n,nbFoc-1)];

vectDist=aux(1,:);
for i=2:n, vectDist=[vectDist aux(i,:)]; end;

card=sum(F');
card(1)=1; % emptyset set to 1
card=repmat(card.^alpha,1,n);

if bal>0
  H=(1-bal)*diag(vectDist.*card/(n*nbFoc))+bal*contraintesMat;
else
  H=diag(vectDist.*card/(n*nbFoc))+contraintesMat;
end
H=sparse(H);



%------------------------ iterations--------------------------------
notfinished=1;
gold=g;

while notfinished

  [masses,lambda,fval]=solqp(H,Aeq,beq,f,x0,0);
  
  x0=masses;

  % reshape m
  m=reshape(masses,nbFoc,n);
  m=m(2:nbFoc,:)';

  % Calculation of centers 
  g=setCentersECM(x,m,F,Smeans,alpha,beta);
  [D,S,Smeans]=setDistances(x,F,g,m,alpha,distance);

  % f matrix
  aux=matContraintes-eye(n);
  contraintesML = max(aux,zeros(n));
  nbContParObjet = sum(contraintesML,2);
  fvide=kron(nbContParObjet,[1; zeros(nbFoc-1,1)]);
  %f=ones(n*nbFoc,1).*x0+fvide;
  f=fvide;

  % H matrix
  aux=D*[zeros(nbFoc-1,1) eye(nbFoc-1)]+[ones(n,1)*rho2, zeros(n,nbFoc-1)];
  vectDist=aux(1,:);
  for i=2:n, vectDist=[vectDist aux(i,:)]; end;

  card=sum(F');
  card(1)=1;
  card=repmat(card.^alpha,1,n);

  H=(1-bal)*diag(vectDist.*card/(n*nbFoc))+bal*contraintesMat;
  H=sparse(H);
  J = masses'*H*masses+bal;

  
  notfinished = sum(sum((abs(g-gold)>1e-3)));
  gold = g ;
end

%--------------------------- end of iterations ----------------------------
m = [abs(ones(n,1)-sum(m,2)) abs(m)];

% BetP calcul
BetP=zeros(n,K);

cardinals = [0 1];
for i = 2:K
  cardinals = [cardinals cardinals+1];
end
cardinals(1) = 1; % to avoid division by 0
aux = m./repmat(cardinals,n,1);

for i = 1:K
  ind = [2^(i-1) + 1 :2^i];
  if i < K
    for j = 1:K-i
      ind = [ind ind+2^(i+j-1)]; 
    end
  end
  % ind = list of indices of supersets of atom i
  BetP(:,i) = sum(aux(:,ind),2);
end


