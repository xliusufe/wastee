#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "_WAST_HEAD.h"


//-------------- Exponential mixture --------------------------------
double Exp_mixture_true(double *y, double *Tn_star, int n, int M, double lambda, double alpha, double *G){
	int i, j, k;
	double tmp, tmp1, tmp2, omega;
	double Tn = 0.0, alpha_inv 	= 1.0/alpha, C1, C2;
	double *ty;


	ty  	= (double*)malloc(sizeof(double)*n);
	for (i = 0; i < n; i++){
		tmp 	= 1.0/(y[i]+lambda);
		ty[i] 	= (y[i] + lambda +1)*exp(y[i]*alpha-y[i])*tmp*tmp;
	}
	for (k = 0; k < M; k++){
		Tn_star[k] = 0.0;
	}
	C1 		= lambda*alpha_inv;
	C2 		= C1*alpha_inv;
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			tmp		= y[i]+y[j]+lambda;
			tmp1 	= tmp + 1;
			tmp2 	= 1.0/tmp;
			omega 	= 1 - C1*(ty[i]+ty[j]) + C2*(tmp1*tmp1+1)*exp((y[i]+y[j])*(alpha-1))*tmp2*tmp2*tmp2;
			for(k=0; k< M; k++){
				Tn_star[k] += omega*G[k*n+i]*G[k*n+j];
			}
			Tn += omega;
		}
	}
	Tn *= 2.0/n;
	for (k = 0; k < M; k++){
		Tn_star[k] *= 2.0/n;
	}

	free(ty);

	return Tn;
}

void Exp_mixture_bstr_true(double *y, double *Tn_star, int n, int M, double lambda, double alpha){
	int i, j, k;
	double tmp, tmp1, tmp2, Tns;
	double alpha_inv, C1, C2;
	double *ty;

	ty  	= (double*)malloc(sizeof(double)*n);

	alpha_inv = 1.0/alpha;
	C1 		= lambda*alpha_inv;
	C2 		= C1*alpha_inv;
	for(k=0; k< M; k++){
		Tns = 0.0;

		for (i = 0; i < n; i++){
			tmp 	= 1.0/(y[k*n+i]+lambda);
			ty[i] 	= (y[k*n+i] + lambda +1)*exp(y[k*n+i]*alpha-y[k*n+i])*tmp*tmp;
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				tmp		= y[k*n+i]+y[k*n+j]+lambda;
				tmp1 	= tmp + 1;
				tmp2 	= 1.0/tmp;
				Tns 	+= - C1*(ty[i]+ty[j]) + C2*(tmp1*tmp1+1)*exp((y[k*n+i]+y[k*n+j])*(alpha-1))*tmp2*tmp2*tmp2;
			}
		}
		Tn_star[k] 	= Tns*2.0/n;
	}


	free(ty);
}

void Exp_mixture_bstr(double *y, double *Tn_star, int n, int M, double lambda, double xlam){
	int i, j, k;
	double tmp, tmp1, tmp2, alpha, Tns;
	double alpha_inv, C1, C2;
	double *ty;

	ty  	= (double*)malloc(sizeof(double)*n);

	for(k=0; k< M; k++){
		Tns 		= 0.0;
		alpha_inv 	= 0.0;
		for (i = 0; i < n; i++){
			alpha_inv += y[k*n+i];
		}
		alpha_inv *= 1.0/n;
		alpha 	= 1.0/alpha_inv;
		C1 		= lambda*alpha_inv;
		C2 		= C1*alpha_inv;

		for (i = 0; i < n; i++){
			tmp 	= 1.0/(y[k*n+i]+lambda);
			ty[i] 	= (y[k*n+i] + lambda +1)*exp(y[k*n+i]*alpha-y[k*n+i])*tmp*tmp;
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				tmp		= y[k*n+i]+y[k*n+j]+lambda;
				tmp1 	= tmp + 1;
				tmp2 	= 1.0/tmp;
				Tns 	+= - C1*(ty[i]+ty[j]) + C2*(tmp1*tmp1+1)*exp((y[k*n+i]+y[k*n+j])*(alpha-1))*tmp2*tmp2*tmp2;
			}
		}
		Tn_star[k] 	= 2.0*Tns/n;
	}


	free(ty);
}

SEXP _Exp_mixture_true(SEXP Y, SEXP Tn_star, SEXP N, SEXP M, SEXP LAMBDA, SEXP ALPHA, SEXP G)
{
	int m = INTEGER(M)[0];
	int n = INTEGER(N)[0];
	SEXP Test;
	PROTECT(Test = allocVector(REALSXP, 1));

	REAL(Test)[0] = Exp_mixture_true(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(ALPHA)[0], REAL(G));

	UNPROTECT(1);
	return Test;
}

SEXP _Exp_mixture_bstr_true(SEXP Y, SEXP N, SEXP M, SEXP LAMBDA, SEXP ALPHA)
{
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Exp_mixture_bstr_true(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(ALPHA)[0]);

	UNPROTECT(1);
	return Tn_star;
}

SEXP _Exp_mixture_bstr(SEXP Y, SEXP N, SEXP M, SEXP LAMBDA, SEXP XLAM)
{
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Exp_mixture_bstr(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(XLAM)[0]);

	UNPROTECT(1);
	return Tn_star;
}

//-------------- Normal mixture --------------------------------
double Normal_mixture_true0(double *y, double *Tn_star, int n, int M, double mu, double sigma2, double alpha, double *G){
	int i, j, k;
	double tmp, tmp1, tmp2, tmp3, omega;
	double Tn = 0.0, C1, C2, C12, C22;
	double *ty;

	C1 	= 1.0/(sigma2+1);
	C2 	= 1.0/(2*sigma2+1);
	C12 = sqrt(C1);
	C22 = sqrt(C2);

	ty  = (double*)malloc(sizeof(double)*n);
	tmp = 0.5*C1;
	for (i = 0; i < n; i++){
		tmp1 	= y[i] - alpha;
		tmp2 	= y[i] - mu;
		ty[i] 	= exp_approx(0.5*tmp1*tmp1 - tmp2*tmp2*tmp, 8);
	}
	for (k = 0; k < M; k++){
		Tn_star[k] = 0.0;
	}
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			tmp 	= y[i]-alpha;
			tmp1 	= y[j]-alpha;
			tmp1 	*= tmp1;
			tmp1 	+= tmp*tmp;
			tmp 	= y[i]-y[j];
			tmp2 	= sigma2*tmp*tmp;
			tmp		= y[i]+y[j]-mu;
			tmp3 	= tmp*tmp;
			omega 	= 1 - C12*(ty[i]+ty[j]) + C22*exp_approx(0.5*tmp1 - 0.5*C2*(tmp2 + tmp3 -2*y[i]*y[j]+mu*mu), 8);
			for(k=0; k< M; k++){
				Tn_star[k] += omega*G[k*n+i]*G[k*n+j];
			}
			Tn += omega;
		}
	}
	Tn *= 2.0/n;
	for (k = 0; k < M; k++){
		Tn_star[k] *= 2.0/n;
	}

	free(ty);

	return Tn;
}

double Normal_mixture_true(double *y, double *Tn_star, int n, int M, double mu, double sigma2, double alpha, double *G){
	int i, j, k;
	double tmp, tmp1, tmp2, tmp3, omega;
	double Tn = 0.0, C1, C2, C12, C22;
	double *ty;

	C1 	= 1.0/(sigma2+1);
	C2 	= 1.0/(2*sigma2+1);
	C12 = sqrt(C1);
	C22 = sqrt(C2);

	ty  = (double*)malloc(sizeof(double)*n);
	tmp = 0.5*C1;
	for (i = 0; i < n; i++){
		tmp1 	= y[i] - alpha;
		tmp2 	= y[i] - mu;
		ty[i] 	= C12*exp(0.5*tmp1*tmp1 - tmp2*tmp2*tmp);
	}
	for (k = 0; k < M; k++){
		Tn_star[k] = 0.0;
	}
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			tmp 	= y[i]-alpha;
			tmp1 	= y[j]-alpha;
			tmp1 	*= tmp1;
			tmp1 	+= tmp*tmp;
			tmp 	= y[i]-y[j];
			tmp2 	= sigma2*tmp*tmp;
			tmp		= y[i]+y[j]-mu;
			tmp3 	= tmp*tmp;
			omega 	= 1 - (ty[i]+ty[j]) + C22*exp(0.5*tmp1 - 0.5*C2*(tmp2 + tmp3 -2*y[i]*y[j]+mu*mu));
			for(k=0; k< M; k++){
				Tn_star[k] += omega*G[k*n+i]*G[k*n+j];
			}
			Tn += omega;
		}
	}
	Tn *= 2.0/n;
	for (k = 0; k < M; k++){
		Tn_star[k] *= 2.0/n;
	}

	free(ty);

	return Tn;
}

void Normal_mixture_bstr(double *y, double *Tn_star, int n, int M, double mu, double sigma2){
	int i, j, k;
	double tmp, tmp1, tmp2, tmp3, omega, alpha, Tns;
	double C1, C2, C12, C22, C3, C4, mu2;
	double *ty;

	ty  	= (double*)malloc(sizeof(double)*n);

	C1 	= 1.0/(sigma2+1);
	C2 	= 1.0/(2*sigma2+1);
	C12 = sqrt(C1);
	C22 = sqrt(C2);
	C3 	= 0.5*C1;
	C4 	= 0.5*C2;
	mu2 = mu*mu;

	for(k=0; k< M; k++){
		Tns 	= 0.0;
		alpha 	= 0.0;
		for (i = 0; i < n; i++){
			alpha += y[k*n+i];
		}
		alpha 	*= 1.0/n;

		for (i = 0; i < n; i++){
			tmp1 	= y[k*n+i] - alpha;
			tmp2 	= y[k*n+i] - mu;
			ty[i] 	= exp(0.5*tmp1*tmp1 - C3*tmp2*tmp2);
		}
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				tmp 	= y[k*n+i]-alpha;
				tmp1 	= y[k*n+j]-alpha;
				tmp1 	*= tmp1;
				tmp1 	+= tmp*tmp;
				tmp 	= y[k*n+i]-y[k*n+j];
				tmp2 	= sigma2*tmp*tmp;
				tmp		= y[k*n+i]+y[k*n+j]-mu;
				tmp3 	= tmp*tmp;
				omega 	= exp(0.5*tmp1 - C4*(tmp2 + tmp3 -2*y[k*n+i]*y[k*n+j]+mu2));
				Tns 	+= - C12*(ty[i]+ty[j]) + C22*omega;
			}
		}
		Tn_star[k] 	= 2.0*Tns/n;
	}

	free(ty);
}

void Normal_mixture_bstr_true(double *y, double *Tn_star, int n, int M, double mu, double sigma2, double alpha){
	int i, j, k;
	double tmp, tmp1, tmp2, tmp3, omega, Tns;
	double C1, C2, C12, C22, C3, C4, mu2;
	double *ty;

	ty  = (double*)malloc(sizeof(double)*n);

	C1 	= 1.0/(sigma2+1);
	C2 	= 1.0/(2*sigma2+1);
	C12 = sqrt(C1);
	C22 = sqrt(C2);
	C3 	= 0.5*C1;
	C4 	= 0.5*C2;
	mu2 = mu*mu;

	for(k=0; k< M; k++){
		Tns 	= 0.0;

		for (i = 0; i < n; i++){
			tmp1 	= y[k*n+i] - alpha;
			tmp2 	= y[k*n+i] - mu;
			ty[i] 	= exp(0.5*tmp1*tmp1 - C3*tmp2*tmp2);
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				tmp 	= y[k*n+i]-alpha;
				tmp1 	= y[k*n+j]-alpha;
				tmp1 	*= tmp1;
				tmp1 	+= tmp*tmp;
				tmp 	= y[k*n+i]-y[k*n+j];
				tmp2 	= sigma2*tmp*tmp;
				tmp		= y[k*n+i]+y[k*n+j]-mu;
				tmp3 	= tmp*tmp;
				omega 	= exp(0.5*tmp1 - C4*(tmp2 + tmp3 -2*y[k*n+i]*y[k*n+j]+mu2));
				Tns 	+= - C12*(ty[i]+ty[j]) + C22*omega;
			}
		}
		Tn_star[k] 	= 2.0*Tns/n;
	}

	free(ty);
}

SEXP _Normal_mixture_true(SEXP Y, SEXP Tn_star, SEXP N, SEXP M, SEXP MU, SEXP SIGMA2, SEXP ALPHA, SEXP G)
{
	int m = INTEGER(M)[0];
	int n = INTEGER(N)[0];
	SEXP Test;
	PROTECT(Test = allocVector(REALSXP, 1));
	if(n>500){
		REAL(Test)[0] = Normal_mixture_true0(REAL(Y), REAL(Tn_star), n, m, REAL(MU)[0], REAL(SIGMA2)[0], REAL(ALPHA)[0], REAL(G));
	}
	else{
		REAL(Test)[0] = Normal_mixture_true(REAL(Y), REAL(Tn_star), n, m, REAL(MU)[0], REAL(SIGMA2)[0], REAL(ALPHA)[0], REAL(G));
	}

	UNPROTECT(1);
	return Test;
}

SEXP _Normal_mixture_bstr(SEXP Y, SEXP N, SEXP M, SEXP MU, SEXP SIGMA2)
{
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Normal_mixture_bstr(REAL(Y), REAL(Tn_star), n, m, REAL(MU)[0], REAL(SIGMA2)[0]);

	UNPROTECT(1);
	return Tn_star;
}

SEXP _Normal_mixture_bstr_true(SEXP Y, SEXP N, SEXP M, SEXP MU, SEXP SIGMA2, SEXP ALPHA)
{
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Normal_mixture_bstr_true(REAL(Y), REAL(Tn_star), n, m, REAL(MU)[0], REAL(SIGMA2)[0], REAL(ALPHA)[0]);

	UNPROTECT(1);
	return Tn_star;
}

//-------------- Poisson mixture --------------------------------
double Poisson_mixture_true(double *y, double *Tn_star, int n, int M, double lambda, double tau, double alpha, double *G){
	int i, j, k;
	double omega, tmp, Tn = 0.0, C1, C2, *ty;

	ty  = (double*)malloc(sizeof(double)*n);

	tmp = exp(alpha)/tgamma(tau);
	C1 	= tmp*pow(lambda/(1.0+lambda),tau);
	C2 	= tmp*exp(alpha)*pow(lambda/(2.0+lambda),tau);

	tmp = alpha+alpha*lambda;
	for (i = 0; i < n; i++){
		ty[i] 	= C1*tgamma(y[i] + tau)/pow(tmp, y[i]);
	}
	for (k = 0; k < M; k++){
		Tn_star[k] = 0.0;
	}

	tmp = 2*alpha+alpha*lambda;
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			omega 	= 1 - (ty[i]+ty[j]) + C2*tgamma(y[i]+y[j]+tau)/pow(tmp, y[i]+y[j]);
			for(k=0; k< M; k++){
				Tn_star[k] += omega*G[k*n+i]*G[k*n+j];
			}
			Tn += omega;
		}
	}
	Tn *= 2.0/n;
	for (k = 0; k < M; k++){
		Tn_star[k] *= 2.0/n;
	}

	free(ty);

	return Tn;
}

void Poisson_mixture_bstr(double *y, double *Tn_star, int n, int M, double lambda, double tau){
	int i, j, k;
	double Tns, tmp, omega, C1, C2, alpha, *ty;

	ty  	= (double*)malloc(sizeof(double)*n);

	for(k=0; k< M; k++){
		Tns 	= 0.0;
		alpha 	= 0.0;
		for (i = 0; i < n; i++){
			alpha += y[k*n+i];
		}
		alpha 	*= 1.0/n;
		tmp = exp(alpha)/tgamma(tau);
		C1 	= tmp*pow(lambda/(1.0+lambda),tau);
		C2 	= tmp*exp(alpha)*pow(lambda/(2.0+lambda),tau);

		tmp = alpha+alpha*lambda;
		for (i = 0; i < n; i++){
			ty[i] 	= C1*tgamma(y[k*n+i] + tau)/pow(tmp, y[k*n+i]);
		}

		tmp = 2*alpha+alpha*lambda;
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				omega 	= tgamma(y[k*n+i]+y[k*n+j]+tau)/pow(tmp, y[k*n+i]+y[k*n+j]);
				Tns 	+= - (ty[i]+ty[j]) + C2*omega;
			}
		}
		Tn_star[k] 	= 2.0*Tns/n;
	}

	free(ty);
}

void Poisson_mixture_bstr_true(double *y, double *Tn_star, int n, int M, double lambda, double tau, double alpha){
	int i, j, k;
	double Tns, tmp, omega, C1, C2, *ty;

	ty  	= (double*)malloc(sizeof(double)*n);

	tmp = exp(alpha)/tgamma(tau);
	C1 	= tmp*pow(lambda/(1.0+lambda),tau);
	C2 	= tmp*exp(alpha)*pow(lambda/(2.0+lambda),tau);
	for(k=0; k< M; k++){
		Tns 	= 0.0;

		tmp = alpha+alpha*lambda;
		for (i = 0; i < n; i++){
			ty[i] 	= C1*tgamma(y[k*n+i] + tau)/pow(tmp, y[k*n+i]);
		}

		tmp = 2*alpha+alpha*lambda;
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				omega 	= tgamma(y[k*n+i]+y[k*n+j]+tau)/pow(tmp, y[k*n+i]+y[k*n+j]);
				Tns 	+= - (ty[i]+ty[j]) + C2*omega;
			}
		}
		Tn_star[k] 	= 2.0*Tns/n;
	}

	free(ty);
}

SEXP _Poisson_mixture_true(SEXP Y, SEXP Tn_star, SEXP N, SEXP M, SEXP LAMBDA, SEXP TAU, SEXP ALPHA, SEXP G){
	int m = INTEGER(M)[0];
	int n = INTEGER(N)[0];
	SEXP Test;
	PROTECT(Test = allocVector(REALSXP, 1));

	REAL(Test)[0] = Poisson_mixture_true(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(TAU)[0], REAL(ALPHA)[0], REAL(G));

	UNPROTECT(1);
	return Test;
}

SEXP _Poisson_mixture_bstr(SEXP Y, SEXP N, SEXP M, SEXP LAMBDA, SEXP TAU){
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Poisson_mixture_bstr(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(TAU)[0]);

	UNPROTECT(1);
	return Tn_star;
}

SEXP _Poisson_mixture_bstr_true(SEXP Y, SEXP N, SEXP M, SEXP LAMBDA, SEXP TAU, SEXP ALPHA)
{
	int n = INTEGER(N)[0];
	int m = INTEGER(M)[0];

	SEXP Tn_star;
	PROTECT(Tn_star = allocVector(REALSXP, m));
	Poisson_mixture_bstr_true(REAL(Y), REAL(Tn_star), n, m, REAL(LAMBDA)[0], REAL(TAU)[0], REAL(ALPHA)[0]);

	UNPROTECT(1);
	return Tn_star;
}

//-------------- Probit regression --------------------------------
void Est_probit_alpha(int *y, double *x, double *alpha, double *residual, int n, int p, int maxIter, double tol){
	int i,j, k,step=0;
	double tmp, tmp1, phix1, phix2, bnorm, bnorm0 = 1.0;
	double *mu, *alpha0, *hess, *dpsi, *psi10;

	mu 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	psi10 	= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alpha0[j] = 1.0;
	}
	while(step<maxIter){
		step++;
		for(j=0; j < p; j++){
			psi10[j] = 0.0;
		}
		AbyB(mu, x, alpha0, n, p, 1);
		for(i=0; i<n; i++){
			tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
			tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
			phix1 	= tmp/tmp1;
			phix2 	= tmp/(1-tmp1);
			tmp 	= y[i]?phix1:(-phix2);
			tmp1 	= y[i]?(-phix1*(mu[i]+phix1)):(phix2*(mu[i] - phix2));
			for(j=0; j < p; j++){
				dpsi[j*n+i] = tmp1*x[j*n+i];
				psi10[j] 	+= tmp*x[j*n+i];
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		MatrixInvSymmetric(hess,p);
    	AbyB(alpha, hess, psi10, p, p, 1);

		bnorm = 0.0;
		for(j=0; j < p; j++){
			alpha[j]	= alpha0[j] - alpha[j];
			bnorm 		+= alpha[j]*alpha[j];
		}

		if(sqrt(bnorm/bnorm0)<tol){
			break;
		}
		else{
			bnorm0 = bnorm;
			for(j=0; j < p; j++){
				alpha0[j] = alpha[j];
			}
		}
	}

	AbyB(mu, x, alpha, n, p, 1);
	for(i=0; i<n; i++){
		tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
		tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
		phix1 	= tmp/tmp1;
		phix2 	= - tmp/(1-tmp1);
		tmp 	= y[i]?phix1:phix2;
		residual[i] = tmp;
	}

	free(mu);
	free(alpha0);
	free(psi10);
	free(dpsi);
	free(hess);
}

void Est_probit1(int *y, double *x, double *alpha, double *residual, double *hess, int n, int p, int maxIter, double tol){
	int i,j, k,step=0;
	double tmp, tmp1, phix1, phix2, bnorm, bnorm0 = 1.0;
	double *mu, *alpha0, *dpsi, *psi10;

	mu 		= (double*)malloc(sizeof(double)*n);
	alpha0 	= (double*)malloc(sizeof(double)*p);
	psi10 	= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0; j < p; j++){
		// alpha0[j] = 0.5;
		alpha0[j] = alpha[j];
	}
	while(step<maxIter){
		step++;
		AbyB(mu, x, alpha0, n, p, 1);
		for(j=0; j < p; j++){
			psi10[j] = 0.0;
		}
		for(i=0; i<n; i++){
			tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
			tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
			phix1 	= tmp/tmp1;
			phix2 	= tmp/(1-tmp1);
			tmp 	= y[i]?phix1:(-phix2);
			tmp1 	= y[i]?(-phix1*(mu[i]+phix1)):(phix2*(mu[i] - phix2));
			for(j=0; j < p; j++){
				dpsi[j*n+i] = tmp1*x[j*n+i];
				psi10[j] 	+= tmp*x[j*n+i];
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}

		MatrixInvSymmetric(hess,p);
    	AbyB(alpha, hess, psi10, p, p, 1);

		bnorm = 0.0;
		for(j=0; j < p; j++){
			alpha[j]	= alpha0[j] - alpha[j];
			bnorm 		+= alpha[j]*alpha[j];
		}
		if( (bnorm<EPS) || (sqrt(bnorm/bnorm0)<tol) ){
			break;
		}
		else{
			bnorm0 = bnorm;
			for(j=0; j < p; j++){
				alpha0[j] = alpha[j];
			}
		}
	}


	AbyB(mu, x, alpha, n, p, 1);
	for(i=0; i<n; i++){
		tmp 	= exp(-0.5*mu[i]*mu[i])*MPI2;
		tmp1 	= 0.5*erf(SQRT2*mu[i]) +0.5;
		phix1 	= tmp/tmp1;
		phix2 	= - tmp/(1-tmp1);
		tmp 	= y[i]?phix1:phix2;
		residual[i] = tmp;
	}

	free(mu);
	free(alpha0);
	free(psi10);
	free(dpsi);
}

void Probit_single_bstr0(double *Tns, int *y, double *tx, double *x, double *z, double *alphahat, int n, int p, int p2,
		int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g, *yb;
	double omega=0.0, Tn, xij;
	double *ty, *resid, *OMEGA, *C11;

	yb 		= (int*)malloc(sizeof(int)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	ty  	= (double*)malloc(sizeof(double)*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alphahat[j] = 0.0;
	}
	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		Est_probit1(yb, tx, alphahat, resid, C11, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(ty);
	free(C11);

}

void Probit_multiple_bstr0(double *Tns, int *y, double *tx, double *x, double *z, double *alphahat, int n, int p, int p2, int p3,
		int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g, *yb;
	double omega=0.0, Tn, rho, sd, xij;
	double *resid, *stdx, *OMEGA, *C11;

	yb 		= (int*)malloc(sizeof(int)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alphahat[j] = 0.0;
	}
	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		Est_probit1(yb, tx, alphahat, resid, C11, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(stdx);
	free(C11);
}

double Probit_single_bstr(double *Tns, int *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2,
		int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g, *yb;
	double omega=0.0, Tn, Tn0=0.0, xij;
	double *ty, *resid, *OMEGA, *C11;

	yb 		= (int*)malloc(sizeof(int)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	ty  	= (double*)malloc(sizeof(double)*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alphahat[j] = 0.0;
	}
	// Est_probit1(y, tx, alphahat, resid, C11, n, p, maxIter, tol);
	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		Est_probit1(yb, tx, alphahat, resid, C11, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(ty);
	free(C11);
	return 2.0*Tn0/n;

}

double Probit_multiple_bstr(double *Tns, int *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2, int p3,
		int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g, *yb;
	double omega=0.0, Tn, Tn0=0.0, rho, sd, xij;
	double *resid, *stdx, *OMEGA, *C11;

	yb 		= (int*)malloc(sizeof(int)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p*p);

	for(j=0; j < p; j++){
		alphahat[j] = 0.5;
	}
	// Est_probit1(y, tx, alphahat, resid, C11, n, p, maxIter, tol);
	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}
	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		Est_probit1(yb, tx, alphahat, resid, C11, n, p, maxIter, tol);
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(stdx);
	free(C11);
	return 2.0*Tn0/n;
}

SEXP _Probit_bstr0(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M;
	double shape1, shape2, tol;

	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn, rAlpha, list, list_names;
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p3==1){
		Probit_single_bstr0(REAL(rTn), INTEGER(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha),
								n, p1, p2, M, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		Probit_multiple_bstr0(REAL(rTn), INTEGER(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha),
								n, p1, p2, p3, M, isBeta, shape1, shape2, maxIter, tol);
	}
	SET_STRING_ELT(list_names, 	0,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	1,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _Probit_bstr(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESID0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M;
	double shape1, shape2, tol;

	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn, rTn0, rAlpha, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p3==1){
		REAL(rTn0)[0] = Probit_single_bstr(REAL(rTn), INTEGER(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
								n, p1, p2, M, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		REAL(rTn0)[0] = Probit_multiple_bstr(REAL(rTn), INTEGER(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
								n, p1, p2, p3, M, isBeta, shape1, shape2, maxIter, tol);
	}
	SET_STRING_ELT(list_names, 	0,  mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	2,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTn);
	SET_VECTOR_ELT(list, 		2, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _Est_probit(SEXP Y, SEXP tX, SEXP DIMs, SEXP TOL){
	int n, p, maxIter;
	n 		= INTEGER(DIMs)[0];
	p 		= INTEGER(DIMs)[1];
	maxIter = INTEGER(DIMs)[2];

	// Outcome
	SEXP rBeta, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	Est_probit_alpha(INTEGER(Y), REAL(tX), REAL(rBeta), REAL(rResids), n, p, maxIter, REAL(TOL)[0]);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

//-------------- Quantile regression --------------------------------
double EstQuantileR1(double *y, double *residual, double tau, const int n){
	double beta=0.0;
	SampleQuantile(&beta, 1, y, n, &tau);

	for(int i=0;i<n;i++){
		residual[i] = y[i] - beta;
	}

	return beta;
}

void EstQuantileR(double *x, double *y, double *beta, double *residual, double qq, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, q2 = 2*qq-1, yk, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);

	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	SampleQuantile(&yk, 1, y, n, &qq);
	for(j=0;j<p;j++)	beta0[j] 	= 0.0;
	for(i=0;i<n;i++)	residual[i]	= y[i] - yk;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			wk = EPS + fabs(residual[i]);
			yk = y[i] + wk*q2;
			for(j=0;j<p;j++){
				tmp 	= x[j*n+i]/wk;
				qy[j] 	+= tmp*yk;
				dpsi[j*n+i] = tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp = beta[j] - beta0[j];
			bnorm += tmp*tmp;
		}

		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)	beta0[j] = beta[j];
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
				residual[i] = y[i] - tmp;
			}
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++) tmp += x[j*n+i]*beta[j];
		residual[i] = y[i] - tmp;
	}

	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
}

void Quantile_single_bstr0(double *Tns, double *y, double *tx, double *x, double *z, double *alphahat, int n, int p, int p2,
		int M, double tau, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, xij;
	double *ty, *resid, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	ty  	= (double*)malloc(sizeof(double)*n);

	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		EstQR2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for(i=0; i<n; i++){
			resid[i] = (IDEX(resid[i], 0) - tau);
		}
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(ty);
}

void Quantile_multiple_bstr0(double *Tns, double *y, double *tx, double *x, double *z, double *alphahat, int n, int p, int p2, int p3,
		int M, double tau, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, rho, sd, xij;
	double *resid, *stdx, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		EstQR2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for(i=0; i<n; i++){
			resid[i] = (IDEX(resid[i], 0) - tau);
		}
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(stdx);
}

double Quantile_single_bstr(double *Tns, double *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2,
		int M, double tau, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, Tn0=0.0, xij;
	double *ty, *resid, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	ty  	= (double*)malloc(sizeof(double)*n);

	// EstQR2(tx, y, alphahat, resid, tau, n, p, maxIter, tol);
	for(i=0; i<n; i++){
		resid0[i] = (IDEX(resid0[i], 0) - tau);
	}
	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstQR2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for(i=0; i<n; i++){
			resid[i] = (IDEX(resid[i], 0) - tau);
		}
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(ty);
	return 2.0*Tn0/n;
}

double Quantile_multiple_bstr(double *Tns, double *y, double *tx, double *x, double *z, double *resid0, double *alphahat, int n, int p, int p2, int p3,
		int M, double tau, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, Tn0=0.0, rho, sd, xij;
	double *resid, *stdx, *yb, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	resid 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	// EstQR2(tx, y, alphahat, resid, tau, n, p, maxIter, tol);
	for(i=0; i<n; i++){
		resid0[i] = (IDEX(resid0[i], 0) - tau);
	}
	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}

			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*resid0[i]*resid0[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstQR2(tx, yb, alphahat, resid, tau, n, p, maxIter, tol);
		for(i=0; i<n; i++){
			resid[i] = (IDEX(resid[i], 0) - tau);
		}
		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*resid[i]*resid[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(OMEGA);
	free(resid);
	free(stdx);
	return 2.0*Tn0/n;
}

SEXP _Quantile_bstr0(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M;
	double tau, shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	tau 	= REAL(PARAMs)[0];
	shape1 	= REAL(PARAMs)[1];
	shape2 	= REAL(PARAMs)[2];
	tol 	= REAL(PARAMs)[3];

	SEXP rTn, rAlpha, list, list_names;
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p3==1){
		Quantile_single_bstr0(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha),
						n, p1, p2, M, tau, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		Quantile_multiple_bstr0(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(rAlpha),
						n, p1, p2, p3, M, tau, shape1, shape2, maxIter, tol);
	}


	SET_STRING_ELT(list_names, 	0,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	1,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _Quantile_bstr(SEXP Y, SEXP tX, SEXP X, SEXP Z, SEXP RESID0, SEXP DIMs, SEXP PARAMs){
	int n, p1, p2, p3, isBeta, maxIter, M;
	double tau, shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p1 		= INTEGER(DIMs)[1];
	p2 		= INTEGER(DIMs)[2];
	p3 		= INTEGER(DIMs)[3];
	M 		= INTEGER(DIMs)[4];
	isBeta 	= INTEGER(DIMs)[5];
	maxIter = INTEGER(DIMs)[6];

	tau 	= REAL(PARAMs)[0];
	shape1 	= REAL(PARAMs)[1];
	shape2 	= REAL(PARAMs)[2];
	tol 	= REAL(PARAMs)[3];

	SEXP rTn, rTn0, rAlpha, list, list_names;
	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p1));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p3==1){
		REAL(rTn0)[0] = Quantile_single_bstr(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
						n, p1, p2, M, tau, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		REAL(rTn0)[0] = Quantile_multiple_bstr(REAL(rTn), REAL(Y), REAL(tX), REAL(X), REAL(Z), REAL(RESID0), REAL(rAlpha),
						n, p1, p2, p3, M, tau, shape1, shape2, maxIter, tol);
	}


	SET_STRING_ELT(list_names, 	0,  mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,  mkChar("Tn"));
	SET_STRING_ELT(list_names, 	2,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTn);
	SET_VECTOR_ELT(list, 		2, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_QR(SEXP Y_, SEXP X_, SEXP PARA_INT_, SEXP PARA_DOUBLE_)
{
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	double *para1 	= REAL(PARA_DOUBLE_);
	double tau   	= para1[0];
	double eps    	= para1[1];

	// Outcome
	SEXP rBeta, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p==1)
		REAL(rBeta)[0]	= EstQuantileR1(REAL(Y_), REAL(rResids), tau, n);
	else
		EstQuantileR(REAL(X_), REAL(Y_), REAL(rBeta), REAL(rResids), tau, n, p, maxstep, eps);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}


//-------------- Double Robust --------------------------------
void EstLogisticR(double *beta, double *residual, double *x, double *y, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *hess, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	hess 	= (double*)malloc(sizeof(double)*p*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*wk*(1.0-wk);
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		wk 	= expx/(1.0+expx);
		yk 	= y[i] - wk;
		residual[i] = yk;
	}

	free(beta0);
	free(hess);
	free(dpsi);
	free(qy);
}

void EstLinearR(double *beta, double *residual, const double *x, const double *y, int n, int p){
	int i,j,k;
	double tmp, *hess, *xy;
	xy 		= (double*)malloc(sizeof(double)*p);
	hess	= (double*)malloc(sizeof(double)*p*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

	MatrixInvSymmetric(hess,p);
	AbyB(beta, hess, xy, p, p, 1);

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0; j < p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(hess);
	free(xy);
}

void EstLogisticR1(double *x, double *y, double *beta, double *residual, double *hess, int n, int p, int maxstep, double eps){
	int i,j,k, step=0;
	double *beta0, *qy, *dpsi;
	double tmp, bnorm, yk, expx, wk;

	beta0 	= (double*)malloc(sizeof(double)*p);
	qy 		= (double*)malloc(sizeof(double)*p);
	dpsi 	= (double*)malloc(sizeof(double)*n*p);

	for(j=0;j<p;j++)	beta0[j] 	= 0.0;

	while (step < maxstep){
		step++;

		for(j=0;j<p;j++) qy[j] = 0.0;
		for(i=0;i<n;i++){
			tmp = 0.0;
			for(j=0;j<p;j++){
				tmp += x[j*n+i]*beta0[j];
			}
			expx = exp(tmp);
			wk 	= expx/(1.0+expx);
			yk 	= wk - y[i];
			tmp = wk*(1.0-wk);
			for(j=0;j<p;j++){
				qy[j] 	+= x[j*n+i]*yk;
				dpsi[j*n+i] = x[j*n+i]*tmp;
			}
		}
		for(j=0; j < p; j++){
			for(k=0; k < p; k++){
				tmp = 0.0;
				for(i=0; i<n; i++){
					tmp += x[j*n+i]*dpsi[k*n+i];
				}
				hess[j*p+k] = tmp;
			}
		}
		MatrixInvSymmetric(hess,p);
    	AbyB(beta, hess, qy, p, p, 1);

		bnorm = 0.0;
		for(j=0;j<p;j++){
			tmp 	=  beta[j];
			bnorm 	+= tmp*tmp;
			beta[j] = beta0[j] - tmp;
		}
		if(sqrt(bnorm)<eps){
			break;
		}
		else{
			for(j=0;j<p;j++)
				beta0[j] = beta[j];
		}
	}

	for(i=0;i<n;i++){
		tmp = 0.0;
		for(j=0;j<p;j++){
			tmp += x[j*n+i]*beta[j];
		}
		expx = exp(tmp);
		wk 	= expx/(1.0+expx);
		yk 	= y[i] - wk;
		residual[i] = yk;
	}

	free(beta0);
	free(qy);
	free(dpsi);
}

void EstLinearR1(double *residual, const double *x, const double *y, double *hess, int n, int p){
    int i,j,k;
    double tmp, *xy, *beta;
    xy 		= (double*)malloc(sizeof(double)*p);
	beta 	= (double*)malloc(sizeof(double)*p);

	for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*n+i]*y[i];
		}
	}
	for(j=0; j < p; j++){
		for(k=0; k < p; k++){
			tmp = 0.0;
			for(i=0; i<n; i++){
				tmp += x[j*n+i]*x[k*n+i];
			}
			hess[j*p+k] = tmp;
		}
	}

    MatrixInvSymmetric(hess,p);
    AbyB(beta, hess, xy, p, p, 1);

	for(i=0; i<n; i++){
		tmp = 0.0;
		for(j=0; j<p; j++){
			tmp += x[j*n+i]*beta[j];
		}
		residual[i] = y[i] - tmp;
	}

	free(beta);
    free(xy);
}

void DoubleRobust_single_bstr0(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *alphahat,
		int n, int p11, int p12, int p2, int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, xij;
	double *resid1, *resid2, *psi, *ty, *C11, *C22, *yb, *ab, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	ab 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	ty 		= (double*)malloc(sizeof(double)*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
			ab[i] = a[g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);
		EstLogisticR1(tx2, ab, alphahat, resid2, C22, n, p12, maxIter, tol);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(ab);
	free(psi);
	free(ty);
	free(resid1);
	free(resid2);
	free(C11);
	free(C22);
	free(OMEGA);
}

void DoubleRobust_multiple_bstr0(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *alphahat,
		int n, int p11, int p12, int p2, int p3, int M, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, rho, sd, xij;
	double *resid1, *resid2, *psi, *stdx, *yb, *ab, *OMEGA, *C11, *C22;

	yb 		= (double*)malloc(sizeof(double)*n);
	ab 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}
			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
			ab[i] = a[g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);
		EstLogisticR1(tx2, ab, alphahat, resid2, C22, n, p12, maxIter, tol);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(ab);
	free(psi);
	free(stdx);
	free(resid1);
	free(resid2);
	free(OMEGA);
	free(C11);
	free(C22);
}

void DoubleRobust_single_bstr1(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *alphahat,
		int n, int p11, int p12, int p2, int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, xij;
	double *resid1, *resid2, *psi, *ty, *C11, *C22, *yb, *ab, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	ab 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	ty 		= (double*)malloc(sizeof(double)*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	EstLogisticR1(tx2, a, alphahat, resid2, C22, n, p12, maxIter, tol);

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(ab);
	free(psi);
	free(ty);
	free(resid1);
	free(resid2);
	free(C11);
	free(C22);
	free(OMEGA);
}

void DoubleRobust_multiple_bstr1(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *alphahat,
		int n, int p11, int p12, int p2, int p3, int M, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, rho, sd, xij;
	double *resid1, *resid2, *psi, *stdx, *yb, *ab, *OMEGA, *C11, *C22;

	yb 		= (double*)malloc(sizeof(double)*n);
	ab 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}
			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
		}
	}

	EstLogisticR1(tx2, a, alphahat, resid2, C22, n, p12, maxIter, tol);
	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(ab);
	free(psi);
	free(stdx);
	free(resid1);
	free(resid2);
	free(OMEGA);
	free(C11);
	free(C22);
}

double DoubleRobust_single_bstr(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *resid01, double *alphahat,
		int n, int p11, int p12, int p2, int M, int isBeta, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, xij, Tn0=0.0;
	double *resid1, *resid2, *psi, *ty, *C11, *C22, *yb, *ab, *OMEGA;

	yb 		= (double*)malloc(sizeof(double)*n);
	ab 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	ty 		= (double*)malloc(sizeof(double)*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);

	EstLinearR1(resid1, tx1, y, C11, n, p11);
	EstLogisticR1(tx2, a, alphahat, resid2, C22, n, p12, maxIter, tol);
	for(i=0; i<n; i++){
		psi[i] 	= resid2[i]*resid01[i];
	}

	if(isBeta==1){
		for(i=0; i<n; i++){
			ty[i] = incbeta(z[i], shape1, shape2);
		}
	}
	else if(isBeta==0){
		shape2 *= SQRT2;
		for(i=0; i<n; i++){
			ty[i] = 0.5*erf( (z[i] - shape1)*shape2 ) +0.5;
		}
	}
	else{
		for(i=0; i<n; i++){
			ty[i] = z[i];
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			omega 	= (z[i]<z[j])?ty[i]:ty[j];
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*psi[i]*psi[j];
		}
	}


	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[g*n+i];
			ab[i] = a[g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(ab);
	free(psi);
	free(ty);
	free(resid1);
	free(resid2);
	free(C11);
	free(C22);
	free(OMEGA);
	return 2.0*Tn0/n;
}

double DoubleRobust_multiple_bstr(double *Tns, double *y, double *a, double *tx1, double *tx2, double *x, double *z, double *resid01, double *alphahat,
		int n, int p11, int p12, int p2, int p3, int M, double shape1, double shape2, int maxIter, double tol){
	int i,j,s,g;
	double omega=0.0, Tn, rho, sd, xij, Tn0=0.0;
	double *resid1, *resid2, *psi, *stdx, *yb, *OMEGA, *C11, *C22;

	yb 		= (double*)malloc(sizeof(double)*n);
	psi  	= (double*)malloc(sizeof(double)*n);
	resid1 	= (double*)malloc(sizeof(double)*n);
	resid2 	= (double*)malloc(sizeof(double)*n);
	stdx 	= (double*)malloc(sizeof(double)*n*p3);
	OMEGA 	= (double*)malloc(sizeof(double)*n*n);
	C11 	= (double*)malloc(sizeof(double)*p11*p11);
	C22 	= (double*)malloc(sizeof(double)*p12*p12);

	// EstLinearR1(resid1, tx1, y, C11, n, p11);
	EstLogisticR1(tx2, a, alphahat, resid2, C22, n, p12, maxIter, tol);
	for(i=0; i<n; i++){
		psi[i] 	= resid2[i]*resid01[i];
	}

	for(i=0; i<n; i++){
		sd = 0.0;
		for(j=0; j<p3; j++){
			sd += z[j*n+i]*z[j*n+i];
		}
		sd = 1.0/sqrt(sd);
		for(j=0; j<p3; j++){
			stdx[j*n+i] = z[j*n+i]*sd;
		}
	}

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			xij = 0.0;
			for (s = 0; s < p2; s++) {
				xij += x[s*n+i]*x[s*n+j];
			}
			rho = 0.0;
			for (s = 0; s < p3; s++) {
				rho += stdx[s*n+i]*stdx[s*n+j];
			}
			if(1-rho*rho < MEPS){
				omega = 0.5;
			}
			else{
				omega = 0.25 + atan(rho/sqrt(1-rho*rho))*MPI1;
			}
			OMEGA[i*n+j]	= omega*xij;
			Tn0 	+= OMEGA[i*n+j]*psi[i]*psi[j];
		}
	}

	for(g=0;g<M;g++){
		Tn = 0.0;
		for(i=0; i<n; i++){
			yb[i] = y[n+g*n+i];
		}
		EstLinearR1(resid1, tx1, yb, C11, n, p11);

		for(i=0; i<n; i++){
			psi[i] 	= resid2[i]*resid1[i];
		}

		for (i = 0; i < n-1; i++) {
			for (j = i+1; j < n; j++) {
				Tn 	+= OMEGA[i*n+j]*psi[i]*psi[j];
			}
		}

		Tns[g] = 2.0*Tn/n;
	}

	free(yb);
	free(psi);
	free(stdx);
	free(resid1);
	free(resid2);
	free(OMEGA);
	free(C11);
	free(C22);
	return 2.0*Tn0/n;
}

SEXP _DoubleRobust_bstr0(SEXP Y, SEXP A, SEXP tX1, SEXP tX2, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p11, p12, p2, p3, M, maxIter, isBeta;
	double shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p11		= INTEGER(DIMs)[1];
	p12		= INTEGER(DIMs)[2];
	p2 		= INTEGER(DIMs)[3];
	p3 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	isBeta 	= INTEGER(DIMs)[6];
	maxIter = INTEGER(DIMs)[7];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn, rAlpha, list, list_names;
  	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p12));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p3==1){
		DoubleRobust_single_bstr0(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(rAlpha),
								n, p11, p12, p2, M, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		DoubleRobust_multiple_bstr0(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(rAlpha),
								n, p11, p12, p2, p3, M, shape1, shape2, maxIter, tol);
	}

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn"));
	SET_STRING_ELT(list_names, 	1,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _DoubleRobust_bstr1(SEXP Y, SEXP A, SEXP tX1, SEXP tX2, SEXP X, SEXP Z, SEXP DIMs, SEXP PARAMs){
	int n, p11, p12, p2, p3, M, maxIter, isBeta;
	double shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p11		= INTEGER(DIMs)[1];
	p12		= INTEGER(DIMs)[2];
	p2 		= INTEGER(DIMs)[3];
	p3 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	isBeta 	= INTEGER(DIMs)[6];
	maxIter = INTEGER(DIMs)[7];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn, rAlpha, list, list_names;
  	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p12));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	if(p3==1){
		DoubleRobust_single_bstr1(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(rAlpha),
								n, p11, p12, p2, M, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		DoubleRobust_multiple_bstr1(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(rAlpha),
								n, p11, p12, p2, p3, M, shape1, shape2, maxIter, tol);
	}

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn"));
	SET_STRING_ELT(list_names, 	1,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn);
	SET_VECTOR_ELT(list, 		1, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _DoubleRobust_bstr(SEXP Y, SEXP A, SEXP tX1, SEXP tX2, SEXP X, SEXP Z, SEXP RESID1, SEXP DIMs, SEXP PARAMs){
	int n, p11, p12, p2, p3, M, maxIter, isBeta;
	double shape1, shape2, tol;
	n 		= INTEGER(DIMs)[0];
	p11		= INTEGER(DIMs)[1];
	p12		= INTEGER(DIMs)[2];
	p2 		= INTEGER(DIMs)[3];
	p3 		= INTEGER(DIMs)[4];
	M 		= INTEGER(DIMs)[5];
	isBeta 	= INTEGER(DIMs)[6];
	maxIter = INTEGER(DIMs)[7];

	shape1 	= REAL(PARAMs)[0];
	shape2 	= REAL(PARAMs)[1];
	tol 	= REAL(PARAMs)[2];

	SEXP rTn0, rTn, rAlpha, list, list_names;
  	PROTECT(rTn0 		= allocVector(REALSXP, 	1));
	PROTECT(rTn 		= allocVector(REALSXP, 	M));
	PROTECT(rAlpha 		= allocVector(REALSXP, 	p12));
	PROTECT(list_names 	= allocVector(STRSXP, 	3));
	PROTECT(list 		= allocVector(VECSXP, 	3));

	if(p3==1){
		REAL(rTn0)[0] = DoubleRobust_single_bstr(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(RESID1), REAL(rAlpha),
								n, p11, p12, p2, M, isBeta, shape1, shape2, maxIter, tol);
	}
	else{
		REAL(rTn0)[0] = DoubleRobust_multiple_bstr(REAL(rTn), REAL(Y), REAL(A), REAL(tX1), REAL(tX2), REAL(X), REAL(Z), REAL(RESID1), REAL(rAlpha),
								n, p11, p12, p2, p3, M, shape1, shape2, maxIter, tol);
	}

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("Tn"));
	SET_STRING_ELT(list_names, 	2,  mkChar("coef"));
	SET_VECTOR_ELT(list, 		0, 	rTn0);
	SET_VECTOR_ELT(list, 		1, 	rTn);
	SET_VECTOR_ELT(list, 		2, 	rAlpha);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(5);
	return list;
}

SEXP _EST_LOGISTICR(SEXP X_, SEXP Y_, SEXP PARA_INT_, SEXP TOL){
	// dimensions
	int *para 	= INTEGER(PARA_INT_);
	int n     	= para[0];
	int p     	= para[1];
	int maxstep	= para[2];

	// Outcome
	SEXP rBeta, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));

	EstLogisticR(REAL(rBeta), REAL(rResids), REAL(X_), REAL(Y_), n, p, maxstep, REAL(TOL)[0]);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}

SEXP _EST_LINEAR(SEXP X_, SEXP Y_, SEXP N, SEXP P){
	int p = INTEGER(P)[0];
	int n = INTEGER(N)[0];
	// Outcome
	SEXP rBeta, rResids, list, list_names;
	PROTECT(rBeta 		= allocVector(REALSXP, 	p));
	PROTECT(rResids		= allocVector(REALSXP, 	n));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2));


	EstLinearR(REAL(rBeta), REAL(rResids), REAL(X_), REAL(Y_), INTEGER(N)[0], p);

	SET_STRING_ELT(list_names, 	0,	mkChar("coef"));
	SET_STRING_ELT(list_names, 	1,  mkChar("residuals"));
	SET_VECTOR_ELT(list, 		0, 	rBeta);
	SET_VECTOR_ELT(list, 		1, 	rResids);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(4);
	return list;
}