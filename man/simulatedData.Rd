\name{simulatedData}
\alias{simulatedData}
\docType{data}
\title{
Simulated data from the framework of general estimating equations
}
\description{
Simulated data from the framework of general estimating equations, including model 
\itemize{
\item 'Gaussian mixture' (\code{simulatedData_mixnorm}), 
\item 'Exponential mixture' (\code{simulatedData_mixexp}), 
\item 'Poission mixture' (\code{simulatedData_mixpoiss}), 
\item 'Quantile regression' (\code{simulatedData_quantile}), 
\item 'Probit regression' (\code{simulatedData_probit}), and 
\item 'Semiparamtric models' (\code{simulatedData_semiparam}).
}

}
\usage{
data(simulatedData_probit)
}
\details{
We simulated data generated from general estimating equations, for instance,
    Probit regression models
    \deqn{
    f(\boldsymbol{V}_i) = \Phi(h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta}))^{Y_i}+\Phi(-h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta}))^{1-Y_i},
    }
    where \eqn{\Phi(\cdot)} is the cumulative distribution function of standard normal distribution, 
    and 
    \deqn{
    h(\boldsymbol{V}_i,\boldsymbol{\alpha}, \boldsymbol{\beta}, \boldsymbol{\theta})
    =\boldsymbol{X}_i\boldsymbol{\alpha}+\boldsymbol{Z}_i^T\boldsymbol{\beta}\boldsymbol{1}(\boldsymbol{U}_i^T\boldsymbol{\theta}\geq0).
    }


    \itemize{
        \item Y: the response, an \eqn{n}-vector

        \item X: the baseline variable with dimension \eqn{n\times p} 

        \item Z: the grouping difference variable with dimension \eqn{n\times q} 

        \item U: the grouping variable with dimension \eqn{n\times r} 
    }
}

\references{
Liu, X. (2022). Testing under loss of identifiability. Manuscript.
}

\examples{
  data(simulatedData_probit)
  y <- data_probit$Y[1:5]
  x <- dim(data_probit$X)
  z <- dim(data_probit$Z)
  u <- dim(data_probit$U)

  data(simulatedData_semiparam)
  y <- data_semiparam$Y[1:5]
  x1 <- dim(data_semiparam$X1)
  x2 <- dim(data_semiparam$X2)
  z <- dim(data_semiparam$Z)
  u <- dim(data_semiparam$U)
}

\keyword{datasets}
