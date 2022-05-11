
EstTn_sup_norm <- function(y, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	rtheta 		= runif(K, 0, 10)
	G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
	if(is.null(alpha)|| isFALSE(alpha))	{
		dims 	= c(n, K, M)
		alpha	= mean(y)
		fit 	<- .Call("_SST_Normal_mixture_bstr", 
							as.numeric(y), 
							as.numeric(rtheta), 
							as.numeric(G), 
							as.integer(dims), 
							as.numeric(alpha)
							)
		pvals 	= fit$pvals
		theta 	= fit$theta
	}
	else{
		Ztheta		= outer(rtheta, y, FUN = function(x, y) exp( (y-alpha)^2/2 -(y-x)^2/2 ) - 1  )
		V_t 		= 1.0 / (exp( (rtheta-alpha)^2) - 1)
		teststat 	= max(rowSums(Ztheta)^2 * V_t)
		G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
		teststat_p 	= apply((Ztheta%*%G)^2 * V_t, 2, max)
		pvals 		= mean(teststat_p >= teststat)
	}
	return(pvals)
}

EstTn_Davies_norm <- function(y, mu = 0, sigma2 = 1, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	if(is.null(alpha)|| isFALSE(alpha))	alpha	= mean(y)

	K = n
	rtheta 		= rnorm(K, mean = mu, sd = sqrt(sigma2))
	rtheta 		= sort(rtheta)
	Ztheta		= outer(rtheta, y, FUN = function(x, y) ( exp( (y-alpha)^2/2 -(y-x)^2/2 ) - 1 )/sqrt(exp((alpha-x)^2) - 1) )
	Ztheta_n 	= rowSums(Ztheta)/sqrt(n)
	teststat 	= max(Ztheta_n)
	V_t 		= sum(abs(Ztheta_n[-1] - Ztheta_n[-K]))

	p_value		= pnorm(-teststat) + V_t*exp(-teststat^2/2)/sqrt(8*pi)

	return(p_value)
}

EstTn_wast_norm <- function(y, mu = 0, sigma2 = 1, alpha = NULL, isBstr = 1, M=1000L) {
	n 		= length(y)
	if(is.null(alpha)|| isFALSE(alpha)){
		if(isBstr==1){
			alpha 		= mean(y)
			teststat 	= .Call("_Normal_mixture_bstr", 
								as.numeric(y), 
								as.integer(n), 
								as.integer(1), 
								as.numeric(mu), 
								as.numeric(sigma2)
								)
			G   		= matrix(data = rnorm(M*n, mean = alpha), nrow = n, ncol = M)
			teststat_p	= .Call("_Normal_mixture_bstr", 
								as.numeric(G), 
								as.integer(n), 
								as.integer(M), 
								as.numeric(mu), 
								as.numeric(sigma2)
								)
		}
		else{
			alpha 		= mean(y)
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Normal_mixture_true", 
								as.numeric(y), 
								as.numeric(teststat_p), 
								as.integer(n), 
								as.integer(M), 
								as.numeric(mu), 
								as.numeric(sigma2), 
								as.numeric(alpha), 
								as.numeric(G)
								)	
		}
	}
	else{
		if(isBstr==1){
			teststat 	<- .Call("_Normal_mixture_bstr_true", 
								as.numeric(y), 
								as.integer(n), 
								as.integer(1), 
								as.numeric(mu), 
								as.numeric(sigma2), 
								as.numeric(alpha)
								)
			alpha 		= mean(y)
			G   		= matrix(rnorm(M*n, mean = alpha), nrow = n, ncol = M)
			teststat_p 	<- .Call("_Normal_mixture_bstr_true", 
								as.numeric(G), 
								as.integer(n), 
								as.integer(M), 
								as.numeric(mu), 
								as.numeric(sigma2), 
								as.numeric(alpha)
								)
		}
		else{
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Normal_mixture_true", 
								as.numeric(y), 
								as.numeric(teststat_p), 
								as.integer(n), 
								as.integer(M), 
								as.numeric(mu), 
								as.numeric(sigma2), 
								as.numeric(alpha), 
								as.numeric(G)
								)	
		}
	}

	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

pval_mixnorm <- function(y, method = "wast", M = 1000, K = 1000, mu = 0, sigma2 = 1, alpha = NULL){
	isBstr 	= 1
	if(method=='sst'){
		pvals  	= EstTn_sup_norm(y, alpha = alpha, K = K, M = M)
	}
	else if(method=='wast') {
	   	pvals  	= EstTn_wast_norm(y, mu = mu, sigma2 = sigma2, alpha = alpha, isBstr = isBstr, M = M)
	}
	else if(method=='davies'){
		pvals 	= EstTn_Davies_norm(y, mu = mu, sigma2 = sigma2, alpha = alpha, K = K, M = M)
	}
	else{
		stop("Method must be one of {'wast', 'sst' and 'davies'} !")
	}

	return(pvals)
}