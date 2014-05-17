[X,Y] = getImageNetXY();

randp = randperm(length(Y));

Y = Y(randp);
X = X(randp,:);

f1arr = getF1Mult(X,Y);