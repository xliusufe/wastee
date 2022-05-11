# wast
R package "wastee" for calculating p-value of the test statistic for subgroup detecting in the framework of general estimating equation (EE). In the paper Liu (2022), we propose a novel U-like statistic by taking the weighted average over the nuisance parametric space. The proposed test statistics not only improve power, but also save dramatically computational time. Many common and useful models are considered, including mixture models and models with change point or change plane. We propose a novel U-like test statistic to detect multiple change planes in the framework of EE.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/wastee")

# Usage

   - [x] [wastee-manual.pdf](https://github.com/xliusufe/wastee/blob/master/inst/wastee-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(wastee)

    data(simulatedData_mixexp)
    pvals <- pval_mixexp(y = data_mixexp, method = "wast")
    pvals

    data(simulatedData_mixnorm)
    pvals <- pval_mixnorm(y = data_mixnorm, method = "wast")
    pvals

    data(simulatedData_mixpoiss)
    pvals <- pval_mixpoiss(y = data_mixpoiss, method = "wast")
    pvals

    data(simulatedData_quantile)
    pvals <- pval_quantile(data = data_quantile, tau = 0.5, method = "wast")
    pvals

    pvals <- pval_quantile(data = data_quantile, tau = 0.3, method = "wast")
    pvals

    data(simulatedData_probit)
    pvals <- pval_probit(data = data_probit, method = "wast")
    pvals

# References
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Davies, R. B. (1977). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 64(2):247-254.

Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 74(1):33-43.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 70(2):291-311.

Liu, X. (2022). Testing under loss of identifiability. Manuscript.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
