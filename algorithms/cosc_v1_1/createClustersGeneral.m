function [allClusters, cut,cheeger,cutPart1,cutPart2,threshold] =  createClustersGeneral(vmin,W,normalized,threshold_type,criterion,deg,notinsideregion)
% Transforms an eigenvector into a cluster indicator function by thresholding.
% 
% Usage:	[allClusters, cut,cheeger,cutPart1,cutPart2,threshold] 
%			= createClusters(vmin,W,normalized,threshold_type,criterion)
%
% Input:
%   vmin: The eigenvector.
%   W: Weight matrix.
%   normalized: True for Ncut/NCC, false for Rcut/RCC.
%   threshold_type: 0: zero, 1: median, 2: mean, -1: best
%   criterion: thresholding criterion. 1: Ratio/Normalized Cut 
%   2: Ratio/Normalized Cheeger Cut
%
%   (optional) deg: Degrees of vertices as column vector. Default is 
%   sum(W,2) in normalized case and ones(size(W,1),1) in unnormalized
%   case. Will be ignored if normalized=false.
%   (optional) notinsideregion=true: don't consider cheeger cuts where you 
%   threshold inside of a region of same value. default is false.
%
% Output:
%   allClusters: Obtained clustering after thresholding.
%   cut: Value of the Normalized/Ratio Cut.
%   cheeger: Value of the Normalized/Ratio Cheeger Cut.
%   cutpart1,cutpart2: The two components of Ratio/Normalized Cut and 
%   Ratio/Normalized Cheeger Cut.
%   threshold: The threshold used to obtain the partitioning.
%
% (C)2010-11 Thomas Buehler and Matthias Hein
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
	
    % Default values for deg
    if (nargin<6)
        if (normalized)
            deg=sum(W,2);
        else
            deg=ones(size(W,1),1);
        end
    end
    %Make deg a row vector;
	if (size(deg,1)>1) 
		deg=deg';
	end

    if (nargin<7)
        notinsideregion=false;
    end

    assert(isempty(find(diag(W)~=0,1)),'Graph contains self loops. W has to have zero diagonal.');
	
	if threshold_type>=0
            threshold= determineThreshold(threshold_type,vmin);
            %allClusters= computeClusterIndicatorFunction(vmin,threshold);
            allClusters= (vmin>threshold);
		    [cutPart1,cutPart2] = computeCutValue(allClusters,W,normalized); %cutPart1: vmin<threshold, cutPart2: vmin>threshold
            cut=cutPart1+cutPart2;
            cheeger=max(cutPart1,cutPart2);
    else

            [vmin_sorted, index]=sort(vmin);
            W_sorted=W(index,index);

            % sum of all degrees in the cluster minus weights within cluster
            deg2=sum(W_sorted); % this has to be the degree also in unnormalized variant
            tempcuts_threshold=cumsum(deg2) - 2*cumsum(full(sum(triu(W_sorted))));
            
            % divide by volume/size
            if(normalized)
                volumes_threshold=cumsum(deg(index));
                 
                cutparts1_threshold=tempcuts_threshold(1:end-1)./volumes_threshold(1:end-1);
                cutparts1_threshold(isnan(cutparts1_threshold))=0;
                cutparts2_threshold=tempcuts_threshold(1:end-1)./(volumes_threshold(end)-volumes_threshold(1:end-1));
                cutparts2_threshold(isnan(cutparts2_threshold))=0;
            else
                sizes_threshold=cumsum(ones(1,size(vmin,1)-1));
                cutparts1_threshold=tempcuts_threshold(1:end-1)./sizes_threshold;
                cutparts2_threshold=tempcuts_threshold(1:end-1)./(size(vmin,1)-sizes_threshold);
            end

            % calculate cuts/cheegers
            cuts_threshold=cutparts1_threshold+cutparts2_threshold;
            cheegers_threshold=max(cutparts1_threshold,cutparts2_threshold);

            
            % also thresholds within regions of same value
            if (~notinsideregion)
            
                % find best cut/cheeger
                if(criterion==1)
                    [cut,threshold_index]=min(cuts_threshold);
                    cheeger=cheegers_threshold(threshold_index);
                elseif(criterion==3)
                    imbalance=0.05;
                    totalsize=size(vmin,1);
                    maxsize=floor(ceil(totalsize/2)*(1+imbalance));
                    [unbalancedcut,threshold_index]=min(tempcuts_threshold(totalsize-maxsize:maxsize));
                    threshold_index=totalsize-maxsize - 1 + threshold_index;
                    cut=cuts_threshold(threshold_index);
                    cheeger=cheegers_threshold(threshold_index);
                else
                    [cheeger,threshold_index]=min(cheegers_threshold);
                    cut=cuts_threshold(threshold_index);
                end

                % update
                cutPart1=cutparts1_threshold(threshold_index);
                cutPart2=cutparts2_threshold(threshold_index);

                allClusters=zeros(size(vmin,1),1);
                allClusters(index(threshold_index+1:end))=1;

                threshold=vmin_sorted(threshold_index);
            
            else
            
                % don't threshold within regions of same value

                %[vminU,indexU]=unique(vmin_sorted(1:end-1));% unique gives index of last occurence
                [vminU,indexU]=unique(vmin_sorted);% unique gives index of last occurence
                vminU=vminU(1:end-1); indexU=indexU(1:end-1); 

                % find best cut/cheeger
                if(criterion==1)
                    [cut,threshold_index]=min(cuts_threshold(indexU));
                    cheeger=cheegers_threshold(indexU(threshold_index));
                elseif(criterion==3) %here nothing has changed
                    imbalance=0.05;
                    totalsize=size(vmin,1);
                    maxsize=floor(ceil(totalsize/2)*(1+imbalance));
                    [unbalancedcut,threshold_index]=min(tempcuts_threshold(totalsize-maxsize:maxsize));
                    threshold_index=totalsize-maxsize - 1 + threshold_index;
                    cut=cuts_threshold(threshold_index);
                    cheeger=cheegers_threshold(threshold_index);
                else
                    [cheeger,threshold_index]=min(cheegers_threshold(indexU));
                    cut=cuts_threshold(indexU(threshold_index));
                end

                % update
                cutPart1=cutparts1_threshold(indexU(threshold_index));
                cutPart2=cutparts2_threshold(indexU(threshold_index));

                threshold=vmin_sorted(indexU(threshold_index));
                allClusters=vmin>threshold;

            end

    end
    
end



function threshold = determineThreshold(threshold_type,u)
% Select the treshold for the cluster indicator function

    assert(threshold_type==0 || threshold_type==1 || threshold_type==2);
    
    switch threshold_type
        case 0
            threshold=0;
        case 1
            threshold = median(u);
        case 2
            threshold = mean(u);
    end
       
end