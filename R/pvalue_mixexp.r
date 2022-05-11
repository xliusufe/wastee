EstTn_sup_exp <- function(y, lambda = 1, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	rtheta 		= runif(K, 1, 10)
	G   		= matrix(data = rnorm(M*n), nrow = n, ncol = M)
	if(is.null(alpha)|| isFALSE(alpha))	{
		alpha	= 1.0/mean(y)
	}
	dims 	= c(n, K, M)
	fit 	<- .Call("_SST_Exp_mixture_bstr",
						as.numeric(y),
						as.numeric(rtheta),
						as.numeric(G),
						as.integer(dims),
						as.numeric(alpha)
						)
	pvals 	= fit$pvals
	theta 	= fit$theta

	return(pvals)
}

EstTn_Davies_exp <- function(y, lambda = 1, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	K = n
	if(is.null(alpha)|| isFALSE(alpha))	alpha	= 1.0/mean(y)

	rtheta 		= rexp(K, rate = lambda) + 1
	rtheta 		= sort(rtheta)

	Ztheta		= outer(rtheta, y, FUN = function(x, y) (x*exp((alpha-x)*y)/alpha - 1)*sqrt(2*x-alpha)*alpha/(x-alpha) )
	Ztheta_n 	= rowSums(Ztheta)/sqrt(n)
	teststat 	= max(Ztheta_n)
	V_t 		= sum(abs(Ztheta_n[-1] - Ztheta_n[-K]))

	p_value 	= pnorm(-teststat) + V_t*exp(-teststat^2/2)/sqrt(8*pi)

	return(p_value)
}

EstTn_wast_exp <- function(y, lambda = 1, alpha = NULL, isBstr = 1, M=1000L) {
	n 		= length(y)
	if(is.null(alpha)|| isFALSE(alpha)){
		if(isBstr==1){
			xlam 		= mean(1/(lambda+rexp(1000000)))
			alpha 		= 1.0/mean(y)
			G   		= matrix(rexp(M*n, rate = alpha), nrow = n, ncol = M)
			teststat 	= .Call("_Exp_mixture_bstr", 
									as.numeric(y),
									as.integer(n), 
									as.integer(1), 
									as.numeric(lambda), 
									as.numeric(xlam)
									)
			teststat_p 	= .Call("_Exp_mixture_bstr", 
									as.numeric(G), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(xlam)
									)
		}
		else{
			alpha 		= 1.0/mean(y);
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Exp_mixture_true", 
									as.numeric(y), 
									as.numeric(teststat_p), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(alpha), 
									as.numeric(G)
									)
		}
	}
	else{
		if (isBstr==1){
			teststat 	= .Call("_Exp_mixture_bstr_true", 
									as.numeric(y), 
									as.integer(n), 
									as.integer(1), 
									as.numeric(lambda), 
									as.numeric(alpha)
									)
			alpha 		= 1.0/mean(y);
			G   		= matrix(rexp(M*n, rate = alpha), nrow = n, ncol = M)
			teststat_p 	= .Call("_Exp_mixture_bstr_true", 
									as.numeric(G), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(alpha)
									)
		}
		else {
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Exp_mixture_true", 
									as.numeric(y), 
									as.numeric(teststat_p), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(alpha), 
									as.numeric(G)
									)
		}
	}

	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

pval_mixexp <- function(y, method = "wast", M = 1000, K = 1000, lambda = 1, alpha = NULL){
	isBstr 	= 1
	if(method=='sst'){
		pvals  	= EstTn_sup_exp(y, lambda = lambda, alpha = alpha, K = K, M = M)
	}
	else if(method=='wast') {
	   	pvals  	= EstTn_wast_exp(y, lambda = lambda, alpha = alpha, isBstr = isBstr, M = M)
	}
	else if(method=='davies'){
		pvals  	= EstTn_Davies_exp(y, lambda = lambda, alpha = alpha, K = K, M = M)
	}
	else{
		stop("Method must be one of {'wast', 'sst' and 'davies'} !")
	}

	return(pvals)
}