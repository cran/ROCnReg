\name{pooledROC.emp}
\alias{pooledROC.emp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical estimation of the pooled ROC curve.
}
\description{
This function estimates the pooled ROC curve using the empirical estimator proposed by Hsieh and Turnbull (1996).
}
\usage{
pooledROC.emp(marker, group, tag.h, data, 
	p = seq(0, 1, l = 101), B = 1000, 
    method = c("ncoutcome", "coutcome"), pauc = pauccontrol(),
  	parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
\arguments{
 \item{marker}{A character string with the name of the diagnostic test variable.}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{Data frame representing the data and containing all needed variables.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the pooled ROC curve.}
  \item{B}{An integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. The default is 1000.}
  \item{method}{A character string specifying if bootstrap resampling (for the confidence intervals) should be done with or without regard to the disease status (``coutcome'' or ``noutcome''). In both cases, a naive bootstrap is used. By default, the resampling is done conditionally on the disease status.}
 \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the pooled ROC curve should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
 \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates the pooled ROC curve (ROC) defined as
\deqn{ROC(p) = 1 - F_{D}\{F_{\bar{D}}^{-1}(1-p)\},}
where
\deqn{F_{D}(y) = Pr(Y_{D} \leq y),}
\deqn{F_{\bar{D}}(y) = Pr(Y_{\bar{D}} \leq y).}
The method implemented in this function estimates \eqn{F_{D}(\cdot)} and \eqn{F_{\bar{D}}(\cdot)} by means of the empirical dsitributions. More precisely, and letting \eqn{\{y_{\bar{D}i}\}_{i=1}^{n_{\bar{D}}}}  and \eqn{\{y_{Dj}\}_{j=1}^{n_{D}}} be two independent random samples from the nondiseased and diseased populations, respectively, the distribution functions in each group take the form
\deqn{\widehat{F}_{D}(y)=\frac{1}{n_D}\sum_{j=1}^{n_D}I(y_{Dj}\leq y),}
\deqn{\widehat{F}_{\bar{D}}(y)=\frac{1}{n_{\bar{D}}}\sum_{i=1}^{n_{\bar{D}}}I(y_{\bar{D}i}\leq y).}

The area under the curve is \deqn{AUC=\int_{0}^{1}ROC(p)dp} and is estimated empirically by means of the Mann-Whitney U-statistic. With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is \deqn{pAUC_{FPF}(u_1)=\int_0^{u_1} ROC(p)dp,} where again is estimated empirically. The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as \deqn{pAUC_{TPF}(u_2)=\int_{u_2}^{1}ROC_{TNF}(p)dp,} where \eqn{ROC_{TNF}(p)} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as \eqn{ROC_{TNF}(p) = F_{\bar{D}}\{F_{D}^{-1}(1-p)\}.} Again, \eqn{ROC_{TNF}(p)} is estimated empirically. The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{marker}{A list with the diagnostic test outcomes in the healthy (h) and diseased (d) groups.}
\item{missing.ind}{A logical value indicating whether missing values occur.}
\item{p}{Set of false positive fractions (FPF) at which the pooled ROC curve has been estimated.}
\item{ROC}{Estimated pooled ROC curve, and corresponding 95\% pointwise confidence band (if computed).}
\item{AUC}{Estimated pooled AUC, and corresponding 95\% confidence interval (if computed).}
\item{pAUC}{If computed, estimated partial area under the pooled ROC curve along with its 95\% confidence interval (if B greater than zero). Note that the returned values are normalised, so that the maximum value is one (see more on Details).}
}
\references{
Hsieh, F., and Turnbull, B.W. (1996). Nonparametric and semiparametric estimation of the receiver operating characteristic curve, The Annals of Statistics, \bold{24}, 25--40.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AROC.bnp}}, \code{\link{AROC.sp}}, \code{\link{AROC.kernel}}, \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}}, \code{\link{pooledROC.dpm}}, \code{\link{cROC.bnp}}, \code{\link{cROC.sp}} or \code{\link{AROC.kernel}}.
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
m0_emp <- pooledROC.emp(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 10,
method = "coutcome", pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"))

summary(m0_emp)

plot(m0_emp)
}
\dontshow{
m0_emp <- pooledROC.emp(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 10,
method = "coutcome", pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"))

summary(m0_emp)

plot(m0_emp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
