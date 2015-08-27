function [clusters,cuts,cheegers, vmin, W, vertex_weights, cnstr1, cnstr2] = ...
    spec_clustering_lcnstrs(W, vertex_weights, CL, ML, normalized, verbosity)

%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
         
    assert(normalized);
    n = size(W,1);
    cheegers = inf;
    cnstr1 = [];
    cnstr2 = [];
    mcnstr1 = [];
    mcnstr2 = [];
    % Remove redundant cannot link constraints!
    if ~isempty(CL)
        dCL = sort(CL, 2);
        dCL = unique(dCL, 'rows');
        q1 = construct_cnstr_graph(CL, n);
        q2 = construct_cnstr_graph(dCL, n);
        assert( sum(sum(q1-q2)) == 0 );
        cnstr1 = dCL(:,1);
        cnstr2 = dCL(:,2);
    end
    
    if ~isempty(ML)
        mcnstr1 = ML(:,1);
        mcnstr2 = ML(:,2);
    end
    totCnstrs = size(cnstr1, 1) + size(mcnstr1,1);
    A = sparse( totCnstrs+1, n );
    
    i=0; j=0;
    if ~isempty(cnstr1)
        for i=1:length(cnstr1)

            if cnstr1(i) ~= cnstr2(i)   % in case of noise, u might have point x cannot link with point x!
                A( i, cnstr1(i) ) = 1;
                A( i, cnstr2(i) ) = 1;
            end

        end
    end

    if ~isempty(mcnstr1)
        for j=1:length(mcnstr1)

            if mcnstr1(j) ~= mcnstr2(j)
                A( j+i, mcnstr1(j) ) = 1;
                A( j+i, mcnstr2(j) ) = -1;
            end

        end
    end
    
    A(i+j+1,:) = vertex_weights';
        
    if verbosity>=3, display('       Solving linearly constrained eigenproblem'); end
    
    if isempty(CL) && isempty(ML)
        [fmin, vmin] = eigs_lcnstrs( W, vertex_weights, A);
%         [vmin, fmin] =eig_std_Laplacian(W,normalized,vertex_weights);
    else
        [fmin, vmin] = eigs_lcnstrs( W, vertex_weights, A);
    end
    
%     if verbosity >= 3
%         fprintf('       eigenvalue  = %.8g\n', fmin);
%     end
    
    if abs(fmin)  <1e-12
       cuts = inf;
       cheegers = inf;
       clusters = ones(n,1);
    else
        res = sign(vmin);
        result = res;
        d = vertex_weights;
        labels = unique(res);
        %    labels
        if length(labels)==2
            cut = sum(sum(W( res==labels(1), res==labels(2) )));

            %d = diag(D);
            vol1 = sum(d(res == labels(1)));
            vol2 = sum(d(res == labels(2)));
            cut = cut/vol1 + cut/vol2;
            cuts = cut(1,1);
            clusters = res;

        else % some vertices have 0 as the value. 

            ix1 = result==1;
            ix2 = result==-1 | result==0;

            cut1 = sum(sum(W( ix1, ix2 )));
            %d = diag(D);
            vol1 = sum(d(ix1));
            vol2 = sum(d(ix2));
            cut1 = cut1/vol1 + cut1/vol2;
            cut1 = cut1(1,1);
            clusters = res;
            clusters(allClusters==0) = -1;

            ix1 = result==-1;
            ix2 = result==1 | result==0;

            cut2 = sum(sum(W( ix1, ix2 )));
            %d = diag(D);
            vol1 = sum(d(ix1));
            vol2 = sum(d(ix2));
            cut2 = cut2/vol1 + cut2/vol2;
            cut2 = cut2(1,1);

            cuts = min(cut1, cut2);
            if cut2 < cut1
                clusters = res;
                clusters(clusters==0) = 1;
            end

        end

    end
    
end