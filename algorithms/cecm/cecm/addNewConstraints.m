function [matConst] = addNewConstraints(x,y,matConst,nbConst, noise, prop);
% x        : The input matrix nxnbAtt.
% y        : The output matrix nx1.
% matConst : The constraints matrix nxn : 1 values is for ML constraints and -1
%            for CL contraits. 0 is when there is no constraint.
% nbConst  : The number of constraint to add.
% noise    : Percentage of noise to include
% prop     : Increase the transitive closure of the constraints if =1.
% @return The constraint matrix with the add of the new contraints


[n,nbAtt]=size(x);

% check the number of constraint to find
nbConstMax=factorial(n-1);
nbConstActual=length(find(matConst-eye(n)~=0));
if nbConstMax<nbConstActual+nbConst
   nbConst=nbConstMax-nbConstActual;
   if sign(nbConst)==-1
      error('The number of constraint is too high !');
   end
end

i=0;
while i<nbConst
 alreadyUsed=0;
 [indMat1 indMat2]=find(matConst-eye(n)~=0);
 ind1=ceil(rand(1)*n);
 ind2=ceil(rand(1)*n);
 while ind1==ind2
   ind2=ceil(rand(1)*n);
 end

 testInd1=[find(indMat1==ind1); find(indMat2==ind1)];
 if ~isempty(testInd1) % if the new indice 1 is already used
   testInd2=[find(indMat2(testInd1)==ind2); find(indMat1(testInd1)==ind2)];
   if ~isempty(testInd2) % if the new pairwise indice is already used
     alreadyUsed=1;
   end
 end

 if ~alreadyUsed %& y(ind1)==y(ind2) % only ML
   matConst(ind1,ind2)=(y(ind1)==y(ind2))*2-1; % constraint setting
   if noise>rand
     matConst(ind1,ind2)=(matConst(ind1,ind2)==-1)*2-1;
   end
   i=i+1;
 end
end

% increasing constraints set
if prop==1
  %% fermeture transitive
  aux1=sign(matConst+matConst'-eye(n));
  aux2=matConst.*sign(matConst);
  aux2=max(aux2,aux2');
  matConst=aux1.*aux2;

  % Must-link propagation
  ML=(matConst>0);
  for i=2:nbConst%n
    ML=ML|ML^i;
  end

  % Cannot-link propagation
  CL=(matConst<0)*1;
  CL=(CL|CL*ML|ML*CL)*1; % dist1
  CL=CL|CL*ML;           % dist2

  matConst=ML-CL;
end
