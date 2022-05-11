\name{exams}
\alias{exams}
\title{
  Examples for subgroup test in the framework of general estimating equation
}
\description{
  Examples of p-value of test statitics for 'Gaussian mixture', 'Exponential mixture', 'Poission mixture', 'Quantile regression', 'Probit regression', and 'Semiparamtric models'.
}

\usage{exams(family = "mixexp", method = "wast", M = 1000, K = 1000, tau = 0.5)
}

\arguments{
  \item{family}{Family for the general estimating equation, including Gaussian mixture ('mixnorm'), Exponential mixture ('mixexp'), Poission mixture ('mixpoiss'), Quantile regression ('quantile'), Probit regression ('probit'), and Semiparamtric models ('semiparam').}
  
  \item{method}{There are two methods, including the proposed 'wast' and 'sst'.}
  
  \item{M}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst'.}

  \item{tau}{The given quantile \eqn{\tau}, a scale in the unit inteval. It is only available if \code{family="quantile"}.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
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

Liu, X. (2022). Testing under loss of identifiability. Manuscript.
}

\examples{ 
  pvals <- exams(family = "mixnorm", method = "wast")
  pvals

  pvals <- exams(family = "mixexp", method = "wast")
  pvals

  pvals <- exams(family = "mixpoiss", method = "wast")
  pvals

  pvals <- exams(family = "quantile", method = "wast", tau = 0.5)
  pvals

  pvals <- exams(family = "probit", method = "wast")
  pvals

  pvals <- exams(family = "semiparam", method = "wast")
  pvals
}