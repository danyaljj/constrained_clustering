function [dML, f] = derive_mls_frm_cls( W, CL )
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    

    dML = cell(0);
    n = size(W,1);
    f = -10*ones(n,1);
    
    if ~isempty(CL)
        u = 1;
        visited = false(n,1);
        
        % Construct constraint graph
        Q = construct_cnstr_graph( CL, n );

        [color, comp] = two_coloring(Q, u, n);
        visited(comp) = true;
        ml1 = comp(color(comp)==1);
        if length(ml1) > 1
            dML = [dML; ml1];
        end
        ml2 = comp(color(comp)~=1);
        if length(ml2) > 1
            dML = [dML; ml2];
        end
        fpartitions(comp) = color(comp);

        
        while sum(visited) < n

            u = find(visited==0, 1);
            [color, comp] = two_coloring(Q, u, n);
            visited(comp) = 1;
            ml1 = comp(color(comp)==1);
            if length(ml1) > 1
                dML = [dML; ml1];
            end
            ml2 = comp(color(comp)~=1);
            if length(ml2) > 1
                dML = [dML; ml2];
            end
            fpartitions(comp) = color(comp);

        end
       
        % Sanity Check
       % assert( sum(fpartitions(CL(:,1)) == fpartitions(CL(:,2))) == 0 );
    end      

end