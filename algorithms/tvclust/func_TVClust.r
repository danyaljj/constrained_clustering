#********************     Variational Inference    *********************
#Created by:    Jeffrey
#Date:          01/17/2013
#Reference:     PRML book by Bishop
#Description:   Variantional Inference for TVClust

#Change:        1. Enforce that V_K = 1. Corrected from VB_DP 01/25/2013
#********************************************************************* 
source("mathFun.r");
 
VB_TVClust = function( x, SM, C, K=6, iterN = 100, isKeepL = 1, alpha0 = 1.2, stopThreshold=0.0005 )
{
	INFINITY = 10^300;
	
	#Begin VB
	#######  Pre Processing
	x = as.matrix(x);
	n = dim(x)[1];
	p = dim(x)[2]; 
	
	meanX = apply(x, 2, mean);
	x = x - rep(meanX, each=n);
	sdX = apply(x, 2, sd);
	x = x / rep(sdX, each=n);

	startT = proc.time();

	#hyper parameters 
	mu0 = rep(0, p); 
	beta0 = 1;
	W0 = 1 * diag(p);
	nu0 = p;
	alpha_p_0 = 1;
	beta_p_0 = 10;
	alpha_q_0 = 10;
	beta_q_0 = 1;

	#initialize parameters
	rQ = matrix(1, n, K)/K; 
	gammaQ = matrix(1, K-1, 2); #para. for stick breaking
	betaQ = rep(beta0, K);   #para for normal wishart dist
	set.seed(1);
	muQ = array(rnorm(p*K), c(p, 1, K));
	WQ = array(0, c(p, p, K));
	for( k in 1 : K )
		WQ[,,k] = W0;
	nuQ = rep(nu0, K);
	alpha_p = 1;
	beta_p  = 10;
	alpha_q = 10;
	beta_q  = 1;
	pImprove = 1;

	#main iterations 
	L = rep(0, iterN);
	nRun = iterN;
	iter = 1;	
	while( iter <= iterN )
	{
		#Find "responsibility" rQ
		for( i in 1 : n )
		{	
			r_i = rep(0, K);
			for( k in 1 : K )
			{
				E_ln_lambda_k = sum( digamma( ( nuQ[k] + 1 - (1:p) ) / 2 ) ) + log( det( WQ[,,k] ) );
				E_sqr = p/betaQ[k] + nuQ[k] * t(x[i,] - muQ[,,k]) %*% WQ[,,k] %*% t(t(x[i,] - muQ[,,k]));
				if( k < K )
					tmp = E_ln_lambda_k / 2 - E_sqr / 2 + digamma( gammaQ[k,1] ) - digamma( sum(gammaQ[k,]) );
				if( k == K )
					tmp = E_ln_lambda_k / 2 - E_sqr / 2;
				if( k > 1 )
					for( j in 1 : (k-1) )
						tmp = tmp + digamma( gammaQ[j,2] ) - digamma( sum(gammaQ[j,]) );
				for( j in 1 : n )
				{
					if( C[i,j] == 1 )
					{	
						tmpj = 0;
						tmpj = tmpj + SM[i,j] * (digamma( alpha_p ) - digamma( alpha_p + beta_p) - digamma( beta_q ) + digamma( alpha_q + beta_q ) ) - (1-SM[i,j]) * (digamma( alpha_q ) - digamma( alpha_q + beta_q ) - digamma( beta_p ) + digamma( alpha_p + beta_p ) );
						tmp = tmp + rQ[j,k] * tmpj; 
					}
				}
				r_i[k] = exp( tmp );
				if( r_i[k] == Inf )
					r_i[k] = INFINITY
			}
			r_i = r_i + 10^(-100); #prevent sum(r_i)=0
			rQ[i,] = r_i / sum(r_i);
		} 
		
		#Update gammaQ
		N = apply(rQ, 2, sum);
		for( k in 1 : (K-1) )
		{
			gammaQ[k,1] = 1 + N[k];
			gammaQ[k,2] = alpha0 + sum(N[(k+1):K]);			
		}
		
		#Update muQ, betaQ, WQ, nuQ
		for( k in 1 : K )
		{
			x_bar_k = apply(x * rQ[,k], 2, sum) / N[k];
			S_k = (t(x) - x_bar_k) %*% ( t(t(x) - x_bar_k) * rQ[,k] ) / N[k];
			
			betaQ[k] = beta0 + N[k];
			muQ[,,k] = (beta0 * mu0 + N[k] * x_bar_k) / betaQ[k];
			tmp = t(t(x_bar_k - mu0)) %*% t(x_bar_k - mu0) * beta0 * N[k] / (beta0 + N[k]);
			WQ[,,k] = solve( tmp + N[k] * S_k + solve(W0) );
			#WQ[,,k] = W0 + tmp + N[k] * S_k; #why not this one?
			nuQ[k] = nu0 + N[k];
		}
		
		#Update alpha_p, beta_p, alpha_q, beta_q  
		tmpRQ = rQ %*% t(rQ); 
		alpha_p = alpha_p_0 + sum( SM * C * tmpRQ );
		beta_p = beta_p_0 + sum( (1-SM) * C * tmpRQ );
		alpha_q = alpha_q_0 + sum( (1-SM) * C * (1-tmpRQ) );
		beta_q = beta_q_0 + sum( SM * C * (1-tmpRQ) );
		
		rQ = rQ + 10^(-40); #prevent log(0)
		rQ = rQ / apply(rQ, 1, sum);
		term75 = sum( rQ * log(rQ) );
		if( term75 > -0.1 )
		{
			print("Optimization converged. Program terminated. @rQ");
			#nRun = iter;
			#iter = iterN+1; #force to quit the loop
		}
		if( (isKeepL == 1) & (iter <= iterN) )
		{   
			#Check convergence by calculating the lower bound (See Eq. 10.71 - 10.77 of Bishop)
			term71 = 0;
			term72 = 0;
			term73 = 0;
			term74 = 0;
			term76 = 0;		
			term77 = 0;
			for( k in 1 : K )
			{
				if( N[k] > 10^(-10) )
				#singularity of 1/N[k]
				{
					E_ln_lambda_k = sum( digamma( ( nuQ[k] + 1 - (1:p) ) / 2 ) ) + log( det( WQ[,,k] ) ); 
					x_bar_k = apply(x * rQ[,k], 2, sum) / N[k];
					S_k = (t(x) - x_bar_k) %*% ( t(t(x) - x_bar_k) * rQ[,k] ) / N[k];
					
					term71 = term71 + .5 * N[k] * ( E_ln_lambda_k - p/betaQ[k] - nuQ[k] * sum( diag( S_k %*% WQ[,,k] ) ) - nuQ[k] * t(x_bar_k - muQ[,,k]) %*% WQ[,,k] %*% t(t(x_bar_k - muQ[,,k])) );
				}
	 
				if( k < K )
				{
					tmp = N[k] * (digamma( gammaQ[k,1] ) - digamma( sum(gammaQ[k,]) ) );	
					for( j in (k+1) : K )
						tmp = tmp + N[j] * ( digamma( gammaQ[k,2] ) - digamma( sum(gammaQ[k,]) ) ); 
					term72 = term72 + tmp;
				}
				
				if( k < K )
					term73 = term73 + (alpha0 - 1) * ( digamma( gammaQ[k,2] ) - digamma( sum(gammaQ[k,]) ) );
				
				tmp2 = ( -p * beta0 / betaQ[k] - beta0 * nuQ[k] * t(muQ[,,k] - mu0) %*% WQ[,,k] %*% t(t(muQ[,,k] - mu0)) ) / 2;
				tmp3 = nuQ[k] * sum( diag( solve( W0 ) %*% WQ[,,k] ) ) / 2;
				term74 = term74 + E_ln_lambda_k * (nu0 - p) / 2 + tmp2 - tmp3;
				
				if( k < K )
					term76 = term76 + (gammaQ[k,1] - 1) * ( digamma( gammaQ[k,1] ) - digamma( sum(gammaQ[k,]) ) ) + (gammaQ[k,2] - 1) * ( digamma( gammaQ[k,2] ) - digamma( sum(gammaQ[k,]) ) ) - lbeta(gammaQ[k,1], gammaQ[k,2]);
				
				term77 = term77 + E_ln_lambda_k / 2 + p/2 * log(betaQ[k]) - H(WQ[,,k], nuQ[k]);
			}  
			
			ElnP = digamma( alpha_p ) - digamma( alpha_p + beta_p);
			ElnMinusP = digamma( beta_p ) - digamma( alpha_p + beta_p);
			ElnQ = digamma( alpha_q ) - digamma( alpha_q + beta_q);
			ElnMinusQ = digamma( beta_q ) - digamma( alpha_q + beta_q);
			tmpRQ = rQ %*% t(rQ);  
			termE_lnPez = sum( SM * C * tmpRQ ) * ElnP + sum( (1-SM) * C * tmpRQ ) * ElnMinusP + sum( SM * C * (1-tmpRQ) ) * ElnMinusQ + sum( (1-SM) * C * (1-tmpRQ) ) * ElnQ;
				
			termE_ln_pq = (alpha_p_0 - alpha_p) * ElnP + (beta_p_0 - beta_p) * ElnMinusP + (alpha_q_0 - alpha_q) * ElnQ + (beta_q_0 - beta_q) * ElnMinusQ + lbeta(alpha_p, beta_p) + lbeta(alpha_q, beta_q);
			
			L[iter] = term71 + term72 + term73 + term74 - term75 - term76 - term77 + termE_lnPez + termE_ln_pq;
			
			pImprove = 1;
			if( iter > 1 )
				pImprove = (L[iter] - L[iter-1]) / abs(L[iter-1]); 
			print(paste("Iter = ", iter, ": L = ", round(L[iter],2), ", improvement = ", round(pImprove*100,4), "%", sep=""));
		}
		if( isKeepL == 0 )
			print(paste("Iter =", iter) );
		if( pImprove < stopThreshold )
		{
			print("Optimization converged. Program terminated. @L");
			nRun = iter;
			iter = iterN + 1; #force to quit the loop
		}
		iter = iter + 1;
	} 
	endT = proc.time();
	runTime = endT - startT;
	zVB = rep(0, n);
	for( i in 1 : n )
		zVB[i] = which.max(rQ[i,]);
 
	#output text:
	print(paste('Running time:', round(runTime[1], 2), 'seconds'));

	result = list(z = zVB, time = runTime[1], L = L, rQ = rQ, gammaQ = gammaQ, nRun=nRun);
}
	
