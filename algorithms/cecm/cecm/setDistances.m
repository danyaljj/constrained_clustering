function [D,Splot,Smean] = setDistances(x,F,g,m,alpha,distance)
  % Calcul des centres agrandis
  [n,nbAtt]=size(x);
  [nbFoc,K]=size(F);

  beta=2;

  gplus=[];
  for i=2:nbFoc
      fi = F(i,:);
      truc = repmat(fi',1,nbAtt);
      gplus = [gplus;sum(g.*truc)./sum(truc)];
  end

  if distance==0
    Splot=repmat({eye(nbAtt)},1,K); % distance euclidienne
    S=repmat({eye(nbAtt)},1,K);
  else % Mahalanobis

    S=[];
    Sigma=[];
    Splot=[]; % matrice de var-cov utilisé pour les plots

    ind=find(sum(F,2)==1);
    for i=1:length(ind) % i=2:nbFoc
      denomSplot=0;
      indi=ind(i);

      Sigmai=zeros(nbAtt,nbAtt);
      for k=1:n % k : le numero de l'individu
	omegai=repmat(F(indi,:),nbFoc,1);
	indAj=find(sum(and(omegai,F),2)>0);

	for j=1:length(indAj)
	  indj=indAj(j);
	  aux = x(k,:)-gplus(indj-1,:);%(x - ones(n,1)*gplus(indj-1,:));
	  Sigmai=Sigmai+sum(F(indj,:)).^(alpha-1)*m(k,indj-1).^beta.*(aux'*aux);

	  denomSplot=denomSplot+sum(F(indj,:)).^(alpha-1)*m(k,indj-1).^beta; % denominateur utilise pour le calcul de Splot (normalisation)
	end
      end

      try
	warning('off');
	lastwarn('','');
	Si=det(Sigmai).^(1/nbAtt)*inv(Sigmai);
	error(lastwarn);
      catch
	Si=det(Sigmai).^(1/nbAtt)*pinv(Sigmai);% dist GK pour Ionosphere
	disp('Utilisation de pinv');
      end

      Splot = [Splot {Sigmai./denomSplot}];
      S = [S {Si}]; % variance des elements singletons uniquement
    end
  end  

  Smean=[];
  for i=1:nbFoc-1
    aux=zeros(nbAtt,nbAtt);
    for j=1:K
      aux=aux+F(i+1,j)*S{j};
    end
    Smean=[Smean {aux./max(sum(F(i+1,:)),1)}];% variance de tous les elements
  end


  % calculation of distances to centers 
  D=[];
  for j=1:nbFoc-1
      aux = (x - ones(n,1)*gplus(j,:));

      if distance==0
	B = diag(aux*aux'); %dist euclidienne
      else
	B = diag(aux*Smean{j}*aux');% dist GK	
      end
      
      D = [D B];
  end % D comprend les distances entre les individus et les centres de gravite (chaque ligne = un individu, chaque colonne = un element sans l'element vide)
%end
