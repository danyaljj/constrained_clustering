source("func_TVClust.r");

############### Generate Synthetic Data ###############
set.seed(1)

# load 
x <- read.table("x2.txt", sep=" ")
C <- read.table("C2.txt", sep=" ")
SM <- read.table("SM2.txt", sep=" ")
#zVB <- read.table("result.txt", sep=" ")
#pObsLink <- read.table("pObsLink.txt", sep=" ")
#trueLables <- read.table("trueLables.txt", sep=" ")

############### Run Program ############### 
K = 15; #initialize K clusters
iterN = 30; #max iterations
isKeepL = 1; #whether monitor the convergence metric
alpha0 = 1; #a tuning para of Dirichlet process
result = VB_TVClust(x, SM, C, K, iterN, isKeepL, alpha0, stopThreshold= 0.00005);

zVB = result$z; 

write.table(zVB, "result2.txt", sep=" ", row.names = FALSE, col.names = FALSE)



