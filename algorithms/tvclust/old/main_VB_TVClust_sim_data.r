source("func_TVClust.r");

############### Generate Synthetic Data ###############
set.seed(1)

sdx = 2;
clusterSize = 30; #number of obs. per cluster
simuData = numeric(0);
for( i in 1 : 1 )
{
	#print(paste("iter:",i));
	x1 = c(rnorm(clusterSize, 5, sdx), rnorm(clusterSize, 0, sdx), rnorm(clusterSize, -5, sdx));
	x2 = c(rnorm(clusterSize, 5, sdx), rnorm(clusterSize, 0, sdx), rnorm(clusterSize, 5, sdx));
	x = cbind(x1, x2); 
} 
trueLables  = rep(1:3, each=clusterSize );

#### Generate side information
n = dim(x)[1];
p = dim(x)[2];
SM0 = (rep(trueLables, each=n) == rep(trueLables, n))*1;
dim(SM0) = c(n,n);
diag(SM0) = 0;
SMpe = 0.05; #SMpe = p(E=1|H=1)
SMqe = 0.95;  #the meaning of SMqe here is different than the one in DP Gibbs sampling code. Here SMqe = p(E=0|H=0)
SM = matrix(-1,n,n);  

pObsLink = 0.1; #pObsLink is the prob of observing a link
tmpC = sample(c(0,1), n*(n-1)/2, replace = T, prob = c(1-pObsLink, pObsLink));
C = matrix(0,n,n);
l = 1;
for( i in 1 : (n-1) )
	for( j in (i+1) : n )
		{
			C[i,j] = tmpC[l];
			C[j,i] = C[i,j];
			l = l + 1;
		} 	

pos_1 = which( (SM0 == 1) & (C==1) );
pos_0 = which( (SM0 == 0) & (C==1) );
SM[pos_1] = sample(c(0,1), length(pos_1), replace = T, prob = c(1-SMpe, SMpe));
SM[pos_0] = sample(c(0,1), length(pos_0), replace = T, prob = c(SMqe, 1-SMqe));
dim(SM) = c(n,n);
diag(SM) = rep(0,n);
for( i in 1 : (n-1) )
	for( j in (i+1) : n )
		SM[j,i] = SM[i,j]; 

#save data
# write.table(x, "x.txt", sep=" ", row.names = FALSE, col.names = FALSE)
# write.table(SM, "SM.txt", sep=" ", row.names = FALSE, col.names = FALSE)
# x = read.table("x.txt");
############### Run Program ############### 
K = 15; #initialize K clusters
iterN = 30; #max iterations
isKeepL = 1; #whether monitor the convergence metric
alpha0 = 1; #a tuning para of Dirichlet process
result = VB_TVClust(x, SM, C, K, iterN, isKeepL, alpha0, stopThreshold= 0.00005);

zVB = result$z; 

write.table(x, "x.txt", sep=" ", row.names = FALSE, col.names = FALSE)
write.table(C, "C.txt", sep=" ", row.names = FALSE, col.names = FALSE)
write.table(SM, "SM.txt", sep=" ", row.names = FALSE, col.names = FALSE)
write.table(zVB, "result.txt", sep=" ", row.names = FALSE, col.names = FALSE)
write.table(pObsLink, "pObsLink.txt", sep=" ", row.names = FALSE, col.names = FALSE)
#write.table(SMpe, "SMpe.txt", sep=" ", row.names = FALSE, col.names = FALSE)
#write.table(SMqe, "SMqe.txt", sep=" ", row.names = FALSE, col.names = FALSE)
write.table(trueLables, "trueLables.txt", sep=" ", row.names = FALSE, col.names = FALSE)

table(zVB, trueLables)



#Visulization 
par(mfrow=c(2,2))
plot(x[,1], x[,2], main="Data: Normal Mixture", type="n", xlab="x", ylab="y");
for( k in 1 : K )
	points( x[trueLables==k,1], x[trueLables==k,2], col=k );
	
plot(x[,1], x[,2], main="Cluster Result: VB for DPMM", type="n", xlab="x", ylab="y");
for( k in 1 : K )
	points( x[zVB==k,1], x[zVB==k,2], col=k );

plot(result$L, main = "Convergence", xlab="Iteration", ylab="Lower Bound L", col=2);
grid(10);


