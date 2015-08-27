logB = function(W, v) #Appendix B.79 of PRML 
{
	p = dim(W)[1];
	result = (-v/2) * log( det(W) );
	tmp = (v*p/2) * log(2) + (p*(p-1)/4) * log(pi) +  sum( lgamma( (v+1 - (1:p) )/2 ) );
	result = result - tmp;
}

H = function( W, v )
{
	p = dim(W)[1];
	E = sum( digamma( (v+1-(1:p)) / 2 ) ) + p * log(2) + log( det(W) ); 
	result = -logB(W,v) - (v-p-1)/2 * E + v * p / 2;
}

CalRandIdx = function(y, SM0)
{
	diag(SM0) = 0;
	n = length(y);
	outputSM = (rep(y, each=n) == rep(y, n))*1;
	dim(outputSM) = c(n,n);
	RandIdx = (sum(SM0 == outputSM)) / n / (n-1); 
}

entropy = function(freqs) 
{
    freqs = freqs/sum(freqs)
    H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
    H = H/log(2)
    return(H)
}
 
GetGeodesicDist = function(DistEuclid, epsilon) 
{
	DistEuclid[DistEuclid > epsilon] = INFINITY;
	g = graph.adjacency(DistEuclid, mode="undirected", weighted=T);
	GeoDist = shortest.paths(g);
}


CalRandIdx = function(y, SM0)
{
	diag(SM0) = 0;
	n = length(y);
	outputSM = (rep(y, each=n) == rep(y, n))*1;
	dim(outputSM) = c(n,n);
	RandIdx = (sum(SM0 == outputSM)) / n / (n-1); 
}


CalRandIdx2 = function(x, y)
{
	n = length(x);
	SMx = (rep(x, each=n) == rep(x, n))*1;
	dim(SMx) = c(n,n);
	SMy = (rep(y, each=n) == rep(y, n))*1;
	dim(SMy) = c(n,n);
	RandIdx = (sum(SMx == SMy) - n) / n / (n-1); 
	RandIdx
}

fillSM = function(SM)
{
	if( !isSymmetric(SM) )
		print("Matrix not symmetric");
	n = dim(SM)[1]; 
	for( i in 1 : (n-2) )
	{
		posOne = which( SM[i, (i+1) : n] == 1 ) + i;
		nOne = length(posOne);
		tmp1 = rep( posOne, times = nOne );
		tmp2 = rep( posOne, each = nOne );
		tmpCombine = cbind( tmp1, tmp2 );
		rows = (tmpCombine[,1] < tmpCombine[,2] );
		if( sum(rows) > 0 )
		{
			index = tmpCombine[rows,];
			if( is.vector(index) )
				index = matrix( index, 1, 2 );
			K = dim(index)[1];
			for( k in 1 : K )
				SM[index[k,1], index[k,2]] = 1; 
		}
	}
	for( i in 1 : (n-1) )
		for( j in (i+1) : n )
			SM[j,i] = SM[i,j];
	diag(SM) = rep(0,n);
	SM;
}	
