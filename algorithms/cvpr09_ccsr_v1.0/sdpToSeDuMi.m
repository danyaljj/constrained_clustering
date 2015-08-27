        
function [A, b, c] = sdpToSeDuMi(F0, FI, cc)
% convert the canonical SDP dual formulation:
% (see  Vandenberche and Boyd 1996, SIAM Review)
%  max -Tr(F0 Z)
% s.t. Tr(Fi Z) = cci and Z is positive definite
%
% in which cc = (cc1, cc2, cc3,..) and FI = {F1, F2, F3,...}
% 
% to SeDuMi format (formulated as vector decision variables ):
% min c'x
% s.t. Ax = b and x is positive definite (x is a vector, so SeDuMi
% really means that vec2mat(x) is positive definite)
%
% by feisha@cis.upenn.edu, June, 10, 2004

if nargin < 3
    error('Cannot convert SDP formulation to SeDuMi formulation in sdpToSeDumi!');
end

[m, n] = size(F0);
if m ~= n
    error('F0 matrix must be squared matrix in sdpToSeDumi(F0, FI, b)');
end

p = length(cc);
if p ~= length(FI)
    error('FI matrix cellarray must have the same length as b in sdpToSeDumi(F0,FI,b)');
end

% should check every element in the cell array FI...later..

% x = reshape(Z, n*n, 1);  % optimization variables from matrix to vector

% converting objective function of the canonical SDP
c = reshape(F0', n*n,1);

% converting equality constraints of the canonical SDP
zz= 0;
for idx=1:length(FI)
    zz= zz + nnz(FI{idx});
end
A = spalloc( n*n, p, zz);
for idx = 1:p
    temp = reshape(FI{idx}, n*n,1);
    lst = find(temp~=0);
    A(lst, idx) = temp(lst);
end
% The SeDuMi solver actually expects the transpose of A as in following
% dual problem
% max b'y
% s.t. c - A'y is positive definite
% Therefore, we transpose A
% A = A';

% b doesn't need to be changed
b = cc;
return;