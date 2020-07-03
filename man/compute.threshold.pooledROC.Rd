\name{compute.threshold.pooledROC}
\alias{compute.threshold.pooledROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pooled ROC based threshold values.
}
\description{
This function implements methods for estimating pooled ROC-based threshold values.
}
\usage{
compute.threshold.pooledROC(object, criterion = c("FPF", "YI"), FPF,
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{pooledROC} as produced by \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}}, or \code{\link{pooledROC.dpm}} functions.}
  \item{criterion}{A character string indicating if the threshold value should be computed based on the Youden index (``YI''), or for fixed false positive fractions (``FPF''). }
  \item{FPF}{For \code{criterion = "FPF"}, a numeric vector with the FPF at which to calculate the threshold values. Atomic values are also valid.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates pooled ROC-based threshold values based on two different criteria, namely, the Youden index (YI) and the one that gives rise to a pre-specified FPF.

The Youden Index is defined as
\deqn{YI = \max_{c}\{TPF(c) - FPF(c)\} = \max_{c}\{F_{\bar{D}}(c) - F_{D}(c)\},}
where
\deqn{F_{D}(y) = Pr(Y_{D} \leq y),}
\deqn{F_{\bar{D}}(y) = Pr(Y_{\bar{D}} \leq y).}
The value \eqn{c^{*}} that achieves the maximum is called the optimal YI threshold. Regarding the criterion for a fixed FPF, the threshold value is obtained as follows
\deqn{c = F_{\bar{D}}^{-1}(1-FPF).}
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{threshold}{If method = "YI", the estimated (optimal) threshold corresponding to the Youden index (the one that maximises TPF/sensitivity + TNF/specificity) and if method = "FPF", the estimated threshold corresponding to the specified FPF. For the Bayesian approaches (\code{pooledROC.dpm} and \code{pooledROC.BB}), and in both cases, in addition to the posterior mean, the 95\% credible interval is also returned.}
\item{YI}{If method = "YI", the estimated Youden index. For the Bayesian approaches (\code{pooledROC.dpm} and \code{pooledROC.BB}), in addition to the posterior mean, the 95\% credible interval is also returned.}
\item{FPF}{If method = "YI",  the FPF corresponding to the estimated (optimal) threshold (For the Bayesian approaches (\code{pooledROC.dpm} and \code{pooledROC.BB}), in addition to the posterior mean, the 95\% credible interval is also returned). If method = "FPF", the supplied FPF argument.}
\item{TPF}{The TPF/sensitivity corresponding to the estimated (optimal) threshold. For the Bayesian approaches (\code{pooledROC.dpm} and \code{pooledROC.BB}), in addition to the posterior mean, the 95\% credible interval is also returned.}
}
\references{
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
\code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}} or \code{\link{pooledROC.dpm}}.
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
m0_dpm <- pooledROC.dpm(marker = "l_marker1", group = "status",
            tag.h = 0, data = newpsa, standardise = TRUE, 
            p = seq(0,1,l=101), compute.WAIC = TRUE, compute.lpml = TRUE, 
            compute.DIC = TRUE, 
            prior.h = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            prior.d = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            mcmc = mcmccontrol(nsave = 400, nburn = 100, nskip = 1))


## Threshold values based on the YI
th_m0_dpm_yi <- compute.threshold.pooledROC(m0_dpm, criterion = "YI")

th_m0_dpm_yi$threshold
th_m0_dpm_yi$YI

### Threshold values for a fixed FPF
th_m0_dpm_fpf <- compute.threshold.pooledROC(m0_dpm, criterion = "FPF", FPF = 0.1)

th_m0_dpm_fpf$threshold
}
\dontshow{
m0_dpm <- pooledROC.dpm(marker = "l_marker1", group = "status",
            tag.h = 0, data = newpsa, standardise = TRUE,  
            p = seq(0,1,l=101), compute.WAIC = FALSE, compute.lpml = FALSE, 
            compute.DIC = FALSE, 
            prior.h = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            prior.d = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            mcmc = mcmccontrol(nsave = 40, nburn = 10, nskip = 1))

## Threshold values based on the YI
th_m0_dpm_yi <- compute.threshold.pooledROC(m0_dpm, criterion = "YI")

th_m0_dpm_yi$threshold
th_m0_dpm_yi$YI

### Threshold values for a fixed FPF
th_m0_dpm_fpf <- compute.threshold.pooledROC(m0_dpm, criterion = "FPF", FPF = 0.1)

th_m0_dpm_fpf$threshold
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
