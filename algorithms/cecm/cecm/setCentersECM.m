function [g]=setCentersECM(x,m,F,Smean,alpha,beta);

[nbFoc,K]=size(F);
[n nbAtt]=size(x);

c = sum(F(2:end,:),2)';  % cardinality of focal sets
indSingleton=find(c==1)+1;

R=[];
B=[];
for l=1:K % pour chaque centre de gravite
  indl=indSingleton(l);

  Rl=[];
  for i=1:n
    Ril=zeros(nbAtt,nbAtt);
    Fl=repmat(F(indl,:),nbFoc,1);
    indAj=find(sum(and(Fl,F),2)==c(indl-1))-1;
    for j=1:length(indAj)
      Ril=Ril+c(indAj(j))^(alpha-1)*m(i,indAj(j))^beta*Smean{indAj(j)};
    end

    Rl=[Rl;Ril];
  end
  R=[R Rl];

  Bl=[];
  for k=1:K
    Bkl=zeros(nbAtt,nbAtt);
    for i=1:n
      indk=indSingleton(k);
      Fl=repmat(sign(F(indl,:)+F(indk,:)),nbFoc,1);
      indAj=find(sum(and(Fl,F),2)==sum(Fl(1,:)))-1;
      for j=1:length(indAj)
	Bkl=Bkl+c(indAj(j))^(alpha-2)*m(i,indAj(j))^beta*Smean{indAj(j)};
      end
    end
    
    Bl=[Bl;Bkl];
  end
  B=[B Bl];
end

X=reshape(x',n*nbAtt,1);
g=B'\R'*X;
g=reshape(g,nbAtt,K)';


