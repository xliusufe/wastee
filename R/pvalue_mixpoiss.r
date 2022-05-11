
EstTn_sup_poiss <- function(y, lambda = 1, tau = 2, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	rtheta 		= runif(K, 0, 10)
	G   		= matrix(data = rnorm(M*n), nrow = n, ncol = M)
	if(is.null(alpha)|| isFALSE(alpha))	{
		dims 	= c(n, K, M)
		alpha	= mean(y)
		fit 	<- .Call("_SST_Poisson_mixture_bstr", 
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
		Ztheta		= outer(rtheta, y, FUN = function(x, y) (x/alpha)^y*exp(alpha-x)-1 )
		V_t 		= 1 /( exp( (alpha-rtheta)^2/alpha ) -1 )
		teststat 	= max(rowSums(Ztheta)^2 * V_t)
		G   		= matrix(data = rnorm(M*n), nrow = n, ncol = M)
		teststat_p 	= apply((Ztheta%*%G)^2 * V_t, 2, max)
		pvals 		= mean(teststat_p >= teststat)
	}

	return(pvals)
}

EstTn_Davies_poiss <- function(y, lambda = 1, tau = 2, alpha = NULL, K=1000L, M=1000L) {
	n 			= length(y)
	K = n
	if(is.null(alpha)|| isFALSE(alpha))	alpha	= mean(y)

	rtheta 		= seq(-10, 10, length.out = K)
	rtheta 		= sort(rtheta)

	Ztheta		= outer(rtheta, y, FUN = function(x, y) ((x/alpha)^y*exp(alpha-x)-1)/sqrt(exp((alpha-x)^2/alpha) - 1) )
	Ztheta_n 	= rowSums(Ztheta)/sqrt(n)
	teststat 	= max(Ztheta_n)
	V_t 		= sum(abs(Ztheta_n[-1] - Ztheta_n[-K]))


	p_value		= pnorm(-teststat) + V_t*exp(-teststat^2/2)/sqrt(8*pi)

	return(p_value)
}

EstTn_wast_poiss <- function(y, lambda = 1, tau = 2, alpha = NULL, isBstr = 1, M=1000L) {
	n 		= length(y)
	if(is.null(alpha)|| isFALSE(alpha)){
		if(isBstr==1){
			alpha		= mean(y)
			G   		= matrix(rpois(M*n, alpha), nrow = n, ncol = M)
			teststat 	<- .Call("_Poisson_mixture_bstr", 
									as.numeric(y), 
									as.integer(n), 
									as.integer(1), 
									as.numeric(lambda), 
									as.numeric(tau)
									)
			teststat_p 	<- .Call("_Poisson_mixture_bstr", 
									as.numeric(G), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(tau)
									)
		}
		else{
			alpha		= mean(y)
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Poisson_mixture_true", 
									as.numeric(y), 
									as.numeric(teststat_p), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(tau), 
									as.numeric(alpha), 
									as.numeric(G)
									)
		}
	}
	else{
		if(isBstr==1){
			teststat 	<- .Call("_Poisson_mixture_bstr_true", 
									as.numeric(y), 
									as.integer(n), 
									as.integer(1), 
									as.numeric(lambda), 
									as.numeric(tau), 
									as.numeric(alpha)
									)
			alpha		= mean(y)
			G   		= matrix(rpois(M*n, alpha), nrow = n, ncol = M)
			teststat_p 	<- .Call("_Poisson_mixture_bstr_true", 
									as.numeric(G), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(tau), 
									as.numeric(alpha)
									)
		}
		else{
			G   		= matrix(rnorm(M*n), nrow = n, ncol = M)
			teststat_p 	= rep(0, M)
			teststat 	<- .Call("_Poisson_mixture_true", 
									as.numeric(y), 
									as.numeric(teststat_p), 
									as.integer(n), 
									as.integer(M), 
									as.numeric(lambda), 
									as.numeric(tau), 
									as.numeric(alpha), 
									as.numeric(G)
									)
		}
	}
	p_value   	= mean(teststat_p >= teststat)

	return(p_value)
}

pval_mixpoiss <- function(y, method = "wast", M = 1000, K = 1000, lambda = 1, tau = 2, alpha = NULL){
	isBstr 	= 1
	if(method=='sst'){
		pvals  	= EstTn_sup_poiss(y, lambda = lambda, tau = tau, alpha = alpha, K = K, M = M)
	}
	else if(method=='wast'){
		pvals  	= EstTn_wast_poiss(y, lambda = lambda, tau = tau, alpha = alpha, isBstr = isBstr, M = M)
	}
	else if(method=='davies'){
		pvals 	= EstTn_Davies_poiss(y, lambda = lambda, tau = tau, alpha = alpha, K = K, M = M)
	}
	else{
		stop("Method must be one of {'wast', 'sst' and 'davies'} !")
	}

	return(pvals)
}