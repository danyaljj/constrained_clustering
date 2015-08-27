function [link] = addConstraints(y,nbConst,noise,option)
% y       : The output matrix nx1.
% nbConst : The number of constraint to add
% noise   : Percentage of noise to include
% option  : Select only one type of constraint
%         : 0 => solely Cannot-Link
%         : 1 => solely Must-Link
%         : 2 => both type of constraints (by default)
% @return The constraints. The two first columns correspond to the
%         pairs of objects, the last column the type of constraint
%         (Cannot-Link=0, Must-Link=1)

n=length(y);

if nargin<3
  noise=0;
end
if nargin<4
  option=2;
end
typeASup=(option-1)*(-1);

link=[];
iConst=nbConst;
while iConst~=0
  link=[link; ceil(rand(iConst,2)*n)];
  [linkinorder idxLink]=unique(sort(link,2),'rows');
  link=link(sort(idxLink),:);
  link(link(:,1)==link(:,2),:)=[]; % suppress links between same points
  typeLink=(y(link(:,1))==y(link(:,2)));
  link(typeLink==typeASup,:)=[]; % suppress link not desired (CL or
                                 % ML or nothing)
  
  iConst=nbConst-size(link,1);
end


if nbConst~=0 
  nbRand=rand(nbConst,1);
  
  link(:,3)=(y(link(:,1))==y(link(:,2)));
  link(nbRand<noise,3)=mod(link(nbRand<noise,3)+1,2);
end
