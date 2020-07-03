\name{compute.threshold.cROC}
\alias{compute.threshold.cROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariate-specific ROC based threshold values.
}
\description{
This function implements methods for estimating covariate-specific ROC-based threshold values.
}
\usage{
compute.threshold.cROC(object, criterion = c("FPF", "YI"), FPF, newdata,
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{cROC} as produced by \code{cROC.bnp()}, \code{cROC.sp()}, or \code{cROC.kernel()}.}
  \item{criterion}{A character string indicating whether the covariate-specific threshold values should be computed based on the Youden index (``YI'') or for fixed set of false positive fractions (``FPF''). }
  \item{FPF}{For \code{criterion = "FPF"}, a numeric vector with the FPF at which to calculate the covariate-specific threshold values. Atomic values are also valid.}
  \item{newdata}{Optional data frame containing the values of the covariates at which the covariate-specific threshold values will be computed. If not supplied, the function \code{\link{cROCData}} is used to build a default dataset.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates covariate-specific ROC-based threshold values based on two different criteria, the Youden index (YI) and the one that gives rise to a pre-specified FPF.

In the conditional case, the Youden index is defined as
\deqn{YI(\mathbf{x}) = \max_{c}|TPF(c|\mathbf{x}) - FPF(c|\mathbf{x})| = \max_{c}|F_{\bar{D}}(c|\mathbf{x}) - F_{D}(c|\mathbf{x})|,}
where
\deqn{F_{D}(y|\mathbf{x}) = Pr(Y_{D} \leq y | \mathbf{X}_{D} = \mathbf{x}),}
\deqn{F_{\bar{D}}(y|\mathbf{x}) = Pr(Y_{\bar{D}} \leq y | \mathbf{X}_{\bar{D}} = \mathbf{x}).}
The value \eqn{c^{*}_{\mathbf{x}}} that achieves the maximum is called the optimal covariate-specific YI threshold. Regarding the criterion for a fixed FPF, the covariate-specific threshold values are obtained as follows
\deqn{c^{*}_{\mathbf{x}} = F_{\bar{D}}^{-1}(1-FPF|\mathbf{x}).}
In both cases, we use the notation \eqn{c^{*}_{\mathbf{x}}} to emphasise that this value depends on covariate \eqn{\mathbf{x}}.
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{newdata}{Data frame containing the values of the covariates at which the covariate-specific thresholds were computed.}
\item{thresholds}{If method = "YI", the estimated covariate-specific (optimal) threshold corresponding to the covariate-specific Youden index (the one that maximises TPF/sensitivity + TNF/specificity) and if method = "FPF", the covariate-specific threshold corresponding to the specified FPF. For the Bayesian approach (\code{cROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
\item{YI}{If method = "YI", the estimated covariate-specific Youden index. For the Bayesian approach (\code{cROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
\item{FPF}{If method = "YI", the FPF corresponding to the estimated (optimal) covariate-specific thresholds (for the Bayesian approach (\code{cROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.). If method = "FPF", the supplied FPF argument.}
\item{TPF}{The covariate-specific TPF/sensitivity corresponding to the estimated covariate-specific (optimal) threshold. For the Bayesian approach (\code{AROC.bnp}), in addition to the posterior mean, the 95\% pointwise credible band is also returned.}
}
\references{
Inacio de Carvalho, V., de Carvalho, M. and Branscum, A. J. (2017). Nonparametric Bayesian Covariate-Adjusted Estimation of the Youden Index. Biometrics, \bold{73}, 1279-1288.

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
\code{\link{cROC.bnp}}, \code{\link{cROC.kernel}} or \code{\link{cROC.sp}}
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
cROC_bnp <- cROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
			  formula.d = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1))

### Threshold values based on the YI
th_cROC_bnp_yi <- compute.threshold.cROC(cROC_bnp, criterion = "YI")

# Plot results
	# Threshold values
	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"est"], 
		type = "l", xlab = "Age", 
		ylab = "log(PSA)", ylim = c(0,3), 
		main = "Threshold values based on the Youden Index")
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"qh"], lty = 2)
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"ql"], lty = 2)

	# Youden Index
	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"est"], 
		type = "l", xlab = "Age", 
		ylab = "log(PSA)", ylim = c(0,1), 
		main = "Threshold values based on the Youden Index")
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"qh"], lty = 2)
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"ql"], lty = 2)

### Threshold values for a fixed FPF
th_cROC_bnp_fpf <- compute.threshold.cROC(cROC_bnp, criterion = "FPF", FPF = 0.1)

# Plot results
	# Threshold values
	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"est"], 
		type = "l", xlab = "Age", 
		ylab = "log(PSA)", ylim = c(0,3), main = "Threshold values for a FPF = 0.1")
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"qh"], lty = 2)
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"ql"], lty = 2)

}
\dontshow{
cROC_bnp <- cROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
			  formula.d = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              mcmc = mcmccontrol(nsave = 10, nburn = 5, nskip = 1))

### Threshold values based on the YI
th_cROC_bnp_yi <- compute.threshold.cROC(cROC_bnp, criterion = "YI")

# Plot results
	# Threshold values
	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"est"], 
		type = "l", xlab = "Age", 
		ylab = "log(PSA)", ylim = c(0,3), 
		main = "Threshold values based on the Youden Index")
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"qh"], lty = 2)
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$thresholds[,"ql"], lty = 2)

	# Youden Index
	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"est"], 
		type = "l", xlab = "Age", 
		ylab = "log(PSA)", ylim = c(0,1), 
		main = "Threshold values based on the Youden Index")
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"qh"], lty = 2)
	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_yi$YI[,"ql"], lty = 2)

### Threshold values for a fixed FPF
#th_cROC_bnp_fpf <- compute.threshold.cROC(cROC_bnp, criterion = "FPF", FPF = 0.1)

## Plot results
#	# Threshold values
#	plot(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"est"], type = "l", xlab = "Age", 
#		ylab = "log(PSA)", ylim = c(0,3), main = "Threshold values for a FPF = 0.1")
#	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"qh"], lty = 2)
#	lines(th_cROC_bnp_yi$newdata$age, th_cROC_bnp_fpf$thresholds[[0.1]][,"ql"], lty = 2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
