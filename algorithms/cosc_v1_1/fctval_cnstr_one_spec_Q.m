function FctVal = fctval_cnstr_one_spec_Q(W, vertex_weights, Q, gamma, fnew, Lovasz)
%
% (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Max Planck Institute for Computer Science, Saarbruecken
% Machine Learning Group, Saarland University, Germany
% http://www.ml.uni-saarland.de
%

%    Q = construct_cnstr_graph(W, CL);
if ~exist('Lovasz','var')
    Lovasz = true;
end
[ix, jx, wval] = find(W);
[qix, qjx, qval] = find(Q);

volQ = sum(sum(Q));
sval = wval.*abs(fnew(ix)-fnew(jx));% + gamma* volQ * (max(fnew) - min(fnew));

if ~Lovasz
    Pfnew = fnew - (fnew'*vertex_weights/sum(vertex_weights));
    FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(vertex_weights'*abs(Pfnew));
else
    n = length(fnew); volV = sum(vertex_weights);
    [fsort,sortind]=sort(fnew);
    sdeg = vertex_weights(sortind);
    if size(sdeg,2)~=1, sdeg = sdeg'; end;
    rcumvols = flipud(cumsum(flipud(sdeg)));
    %rcumvols2 = sum(sdeg) - [0; cumsum(sdeg(1:n-1))];
    %assert( sum(abs(rcumvols-rcumvols2)) == 0 );
    vec = zeros(n,1);
    vec(sortind) = 2*sdeg.*(volV - rcumvols - [rcumvols(2:end); 0])/volV;
    FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(vec'*fnew);  % we cancelled factor 0.5 from numerator and denominator.
end
end
