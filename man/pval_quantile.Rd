\name{pval_quantile}
\alias{pval_quantile}
\title{
  P-value for subgroup test in quantile regression models
}
\description{
  Provide p-value for subgroup test in quantile regression models, including two methods 'wast' and 'sst'.
}

\usage{
  pval_quantile(data, method = "wast", tau = 0.5, M = 1000, K = 2000, 
                              isBeta = FALSE, shape1 = 1, shape2 = 1)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}
  
  \item{method}{There are two methods, including the proposed 'wast' and 'sst'.}

  \item{tau}{The given quantile \eqn{\tau}, a scale in the unit inteval.}

  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst'.}

  \item{isBeta}{A bool value. The weight \eqn{w(\gamma)} is chosen to be Beta distribution if \code{isBeta=TRUE}, which can be used if the grouping difference variable is bounded in \eqn{[0,1]}. Default is \code{FALSE}.}

  \item{shape1}{The first parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{shape2}{The second parameter of Best distribution if \code{isBeta = TRUE}.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Quantile regression models
  \deqn{
  Q_{Y_i}(\tau|\boldsymbol{X}_i,\boldsymbol{Z}_i, \boldsymbol{U}_i) 
  = \boldsymbol{X}_i^T\boldsymbol{\alpha}(\tau)+\boldsymbol{Z}_i^T\boldsymbol{\beta}(\tau)\boldsymbol{ 1}(\boldsymbol{U}_i^T\boldsymbol{\theta}(\tau)\geq0).
  }

The hypothesis test problem is 
\deqn{
  H_0: \boldsymbol{\beta}=\boldsymbol{0}\quad versus\quad H_1: \boldsymbol{\beta}\neq \boldsymbol{0}.
}
}

\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

LEE, S., SEO, M. H. and SHIN, Y. (2011). Testing for Threshold Effects in Regression Models. Journal of
the American Statistical Association 106, 220-231.

Liu, X. (2022). Testing under loss of identifiability. Manuscript.
}

\examples{ 
  data(simulatedData_quantile)
  pvals <- pval_quantile(data = data_quantile, tau = 0.5, method = "wast")
  pvals

  pvals <- pval_quantile(data = data_quantile, tau = 0.3, method = "wast")
  pvals

}