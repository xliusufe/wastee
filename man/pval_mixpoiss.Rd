\name{pval_mixpoiss}
\alias{pval_mixpoiss}
\title{
  P-value for subgroup test in Poission mixture models
}
\description{
  Provide p-value for subgroup test in Poission mixture models, including two methods 'wast' and 'sst'.
}

\usage{
  pval_mixpoiss(y, method = "wast", M = 1000, K = 1000, 
                      lambda = 1, tau = 2, alpha = NULL)
}

\arguments{
  \item{y}{The data.}
  
  \item{method}{There are three methods, including the proposed 'wast', 'sst', and 'davies'. Refer method 'davies' to Davies (1977, 1987)
      which is only avaliable for known true parameter \eqn{\alpha=\alpha_0}.}

  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'davies'.}

  \item{lambda}{The rate parameter \eqn{\lambda} of weight that is a Gamma distribution.}

  \item{tau}{The shape parameter \eqn{\tau} of weight that is a Gamma distribution.}

  \item{alpha}{Mean under the null. Default is \code{NULL}, in which the mean \eqn{\alpha} is eatimeted by the sample mean.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Poisson mixture model
  \deqn{
  f(y) = (1-\beta)\frac{\alpha^{y}\exp(-\alpha)}{y!} + \beta\frac{\theta^{y}\exp(-\theta)}{y!},
  }
  where \eqn{\beta\in[0,1]},  \eqn{\alpha\in\boldsymbol{R}} and \eqn{\theta\in\boldsymbol{R}} are unknown parameters.
}

\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Davies, R. B. (1977). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 64(2):247-254.

Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 74(1):33-43.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

LEE, S., SEO, M. H. and SHIN, Y. (2011). Testing for Threshold Effects in Regression Models. Journal of
the American Statistical Association 106, 220-231.

Liu, X. (2022). Testing under loss of identifiability. Manuscript.
}

\examples{ 
  data(simulatedData_mixpoiss)
  pvals <- pval_mixpoiss(y = data_mixpoiss, method = "wast")
  pvals

  pvals <- pval_mixpoiss(y = data_mixpoiss, method = "sst")
  pvals

  pvals <- pval_mixpoiss(y = data_mixpoiss, method = "davies", alpha = 1)
  pvals

}