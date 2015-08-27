function dist = euclidean(X1, X2)
    L1 = X1  - X2;
    dist = sqrt(L1 * L1');
end