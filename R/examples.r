exams <- function(family = "mixexp", method = "wast", M = 1000, K = 1000, tau = 0.5){
	if(family == "mixexp"){
		data(simulatedData_mixexp)
		pvals   = pval_mixexp(y = data_mixexp, method = method, M=M, K = K)
	}
	else if(family == "mixpoiss"){
		data(simulatedData_mixpoiss)
		pvals   = pval_mixpoiss(y = data_mixpoiss, method = method, M=M, K = K)
	}
	else if(family == "mixnorm"){
		data(simulatedData_mixnorm)
		pvals   = pval_mixnorm(y = data_mixnorm, method = method, M=M, K = K)
	}
	else if(family == "probit"){
		data(simulatedData_probit)
		pvals   = pval_probit(data = data_probit, method = method, M=M, K = K)
	}
	else if(family == "quantile"){
		data(simulatedData_quantile)
		pvals   = pval_quantile(data = data_quantile, method = method, tau = tau, M=M, K = K)
	}
	else if(family == "semiparam"){
		data(simulatedData_semiparam)
		pvals   = pval_semiparam(data = data_semiparam, method = method, M=M, K = K)
	}
	else{
		stop("Family must be one of {Gaussian mixture ('mixnorm'), Exponential mixture ('mixexp'), Poission mixture ('mixpoiss'), Quantile regression ('quantile'), Probit regression ('probit'), and Semiparamtric models ('semiparam')} !")
	}

	return(pvals)
}