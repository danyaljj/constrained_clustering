function fsolns = top_fsolns( W, CL, Y, C )
%% the kth feasible solution, fsolns(:, k) has three values: 0, 1, -10. 
%% where -10 is for unconstrained vertices.
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%    
    n = size(W,1);
    f = -10*ones(n,1);
    Q = sparse(CL(:,1), CL(:,2), 1, n,n) + sparse(CL(:,2), CL(:,1), 1, n,n);
    
    % find connected components of Q
    % for each connected component do the breadth first traversal assigning
    % labels

    [comp,connected,sizes]=connectedComponents(Q);
    
    components = find(sizes>1);
    component_sizes = sizes(components);
    [scs, s_ix] = sort(component_sizes, 'descend');
%    [components, 
    jx = [];
    ff = [];
    counter = 0;
    for i=1:length(components)
       
        ix = find(comp==components(s_ix(i)));
        %ix = find(comp==components(i));
        Qcomp = Q(ix, ix);
        [labels] = bfs_traversal(Qcomp);
        f(ix) = labels;
        
    end
   
   display(counter);
   %display(length(components));
   assert( sum(f(CL(:,1)) == f(CL(:,2))) == 0 );
   
   %f(f==0)=2;
   
   %fsolns = -10*ones(2^k,1);
   C = min(C, length(components));
   fsolns = -10*ones(n,2^C);
   fsolns(:,end) = f;
   
   y = unique(Y);
   f( f==0 ) = y(1);
   f( f==1 ) = y(2);
   
  assert( sum(f(CL(:,1)) == Y(CL(:,1))) == sum(f(CL(:,1)) == Y(CL(:,1))) );
  display( max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) ));

  sats = max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) );
%    for i=1:length(components)
%        
%        fsolns(:,i) = fsolns(:,1);
%        ix = find(comp==components(s_ix(i)));
%        fsolns(ix,i) = ~fsolns(ix,1);
%        
% %        ix = find(comp==components(s_ix(i+1)));
% %        fsolns(ix,i) = ~fsolns(ix,1);
% 
% %         ix = find(comp==components(s_ix(i+2)));
% %        fsolns(ix,i) = ~fsolns(ix,1);
% 
%        f = fsolns(:,i);
%        y = unique(Y);
%        f( f==0 ) = y(1);
%        f( f==1 ) = y(2);
% 
%        assert( sum(f(CL(:,1)) == Y(CL(:,1))) == sum(f(CL(:,1)) == Y(CL(:,1))) );
%        display( max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) ));
%        
%        sats = max(sats, max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) ));
%        
%    end

%    C = 5;

    for i=1:2^C-1
       
   
        b = dec2bin(i);
        b = b=='1';
        l = size(b,2);
        labels(1:C-l) = 0;
        labels(C-l+1:C) = b;


        fsolns(:, i) = fsolns(:,end);
        for k=1:C
    
            ix = find(comp==components(s_ix(k)));
            fsolns(ix,i) = xor( fsolns(ix,1), labels(k) ) ;
            
        end
        counter = counter+1;
       
       f = fsolns(:,i);
       y = unique(Y);
       f( f==0 ) = y(1);
       f( f==1 ) = y(2);

       assert( sum(f(CL(:,1)) == Y(CL(:,1))) == sum(f(CL(:,1)) == Y(CL(:,1))) );
%        display( max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) ));
       
       sats = max(sats, max (sum(f(CL(:,1)) == Y(CL(:,1))), sum(f(CL(:,1)) ~= Y(CL(:,1))) ));

    end
   %display(length(components));
%    counter
%    sats
end