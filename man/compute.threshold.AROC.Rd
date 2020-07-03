\name{compute.threshold.AROC}
\alias{compute.threshold.AROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AROC based threshold values.
}
\description{
This function implements methods for estimating AROC-based threshold values.
}
\usage{
compute.threshold.AROC(object, criterion = c("FPF", "YI"), FPF, newdata,
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{AROC} as produced by \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()}.}
  \item{criterion}{A character string indicating whether the covariate-adjusted threshold values should be computed based on the Youden index (``YI'') or for a fixed set of false positive fractions (``FPF'').}
  \item{FPF}{For criterion = FPF, a numeric vector with the FPF at which to calculate the AROC-based threshold values. Atomic values are also valid.}
  \item{newdata}{Optional data frame containing the values of the covariates at which the AROC-based threshold values will be computed. If not supplied, the function \code{\link{cROCData}} is used to build a default dataset.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates AROC-based threshold values based on two different criteria, namely, the Youden index (YI) and the one that gives rise to a pre-specified FPF. Before proceeding, we would like to mention that when the accuracy of a test is not affected by covariates, this does not necessarily imply that the covariate-specific ROC curve (which in this case is the same for all covariate values) coincides with the pooled ROC curve. It does coincide, however, with the AROC curve. Consequently, in all cases where covariates affect the test, even though they might not affect its discriminatory capacity, inferences based on the pooled ROC curve might be misleading. In such cases the AROC curve should be used instead. This also applies to the selection of (optimal) threshold values, which, as will be seen, might be covariate-specific (i.e., possibly different for different covariate values).

For the AROC curve, the Youden Index is defined as
\deqn{YI = \max_{p}\{AROC(p) - p\},}
The value \eqn{p^{*}} (FPF) that achieves the maximum is then used to calculate the optimal (covariate-specific) YI threshold as follows
\deqn{c^{*}_{\mathbf{x}} = F_{\bar{D}}^{-1}(1-p^{*}|\mathbf{x}),}
where
\deqn{F_{\bar{D}}(y|\mathbf{x}) = Pr(Y_{\bar{D}} \leq y | \mathbf{X}_{\bar{D}} = \mathbf{x}).}
In a similar way, when using the criterion for a fixed FPF, the covariate-specific threshold values are obtained as follows
\deqn{c^{*}_{\mathbf{x}} = F_{\bar{D}}^{-1}(1-FPF|\mathbf{x}).}
In both cases, we use the notation \eqn{c^{*}_{\mathbf{x}}} to emphasise that this value depends on covariate \eqn{\mathbf{x}}.
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{newdata}{Data frame containing the values of the covariates at which the AROC-based thresholds were computed.}
\item{thresholds}{If method = "YI", the estimated AROC-based threshold corresponding to the Youden index, and if method = "FPF", AROC-based threshold corresponding to the specified FPF. For the Bayesian approach (\code{AROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
\item{YI}{If criterion = "YI", the AROC-based Youden index. For the Bayesian approach (\code{AROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
\item{FPF}{If criterion = "YI", the FPF where the Youden index is attained, and if criterion = "FPF", the supplied FPF argument. For the Bayesian approach (\code{AROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
}
\references{
Inacio de Carvalho, V., and Rodriguez-Alvarez, M. X. (2018). Bayesian nonparametric inference for the covariate-adjusted ROC curve. arXiv preprint arXiv:1806.00473.

Rodriguez-Alvarez, M. X., Roca-Pardinas, J., and Cadarso-Suarez, C. (2011). ROC curve and covariates: extending induced methodology to the non-parametric framework. Statistics and Computing, \bold{21}, 483--499.

Rutter, C.M. and Miglioretti, D. L. (2003). Estimating the Accuracy of Psychological Scales Using Longitudinal Data. Biostatistics, \bold{4}, 97--107.

Youden, W. J. (1950). Index for rating diagnostic tests. Cancer, \bold{3}, 32--35.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{AROC.bnp}}, \code{\link{AROC.kernel}} or \code{\link{AROC.sp}}
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              prior = priorcontrol.bnp(m0 = rep(0, 4), 
              S0 = 10*diag(4), nu = 6, Psi = diag(4),
              a = 2, b = 2, aalpha = 2, balpha = 2, L =10),
              mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1))

### Threshold values based on the YI
th_AROC_bnp_yi <- compute.threshold.AROC(AROC_bnp, criterion = "YI")

# Plot results
plot(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"est"], 
	type = "l", xlab = "Age", 
	ylab = "log(PSA)", ylim = c(0,3), 
	main = "Threshold values based on the Youden Index")
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"qh"], lty = 2)
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"ql"], lty = 2)

### Threshold values for a fixed FPF
th_AROC_bnp_fpf <- compute.threshold.AROC(AROC_bnp, criterion = "FPF", FPF = 0.1)

# Plot results
plot(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"est"], 
	type = "l", xlab = "Age", 
	ylab = "log(PSA)", ylim = c(0,3), 
	main = "Threshold values for a FPF = 0.1")
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"qh"], lty = 2)
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"ql"], lty = 2)

}
\dontshow{
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              prior = priorcontrol.bnp(m0 = rep(0, 4), S0 = 10*diag(4), nu = 6, Psi = diag(4),
              a = 2, b = 2, aalpha = 2, balpha = 2, L =10),
              mcmc = mcmccontrol(nsave = 20, nburn = 5, nskip = 1))

### Threshold values based on the YI
th_AROC_bnp_yi <- compute.threshold.AROC(AROC_bnp, criterion = "YI")

# Plot results
plot(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"est"], 
	type = "l", xlab = "Age", 
	ylab = "log(PSA)", ylim = c(0,3), 
	main = "Threshold values based on the Youden Index")
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"qh"], lty = 2)
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_yi$thresholds[,"ql"], lty = 2)

### Threshold values for a fixed FPF
th_AROC_bnp_fpf <- compute.threshold.AROC(AROC_bnp, criterion = "FPF", FPF = 0.1)

# Plot results
plot(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"est"], 
	type = "l", xlab = "Age", 
	ylab = "log(PSA)", ylim = c(0,3), 
	main = "Threshold values for a FPF = 0.1")
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"qh"], lty = 2)
lines(th_AROC_bnp_yi$newdata$age, th_AROC_bnp_fpf$thresholds[[0.1]][,"ql"], lty = 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
