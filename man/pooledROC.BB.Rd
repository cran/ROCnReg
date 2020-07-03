\name{pooledROC.BB}
\alias{pooledROC.BB}
\title{
Bayesian bootstrap estimation of the pooled ROC curve.
}
\description{
This function estimates the pooled ROC curve using the Bayesian bootstrap estimator proposed by Gu et al. (2008).
}
\usage{
pooledROC.BB(marker, group, tag.h, data, 
	p = seq(0, 1, l = 101), B = 5000, pauc = pauccontrol(),
  	parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{A character string with the name of the diagnostic test variable.}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{Data frame representing the data and containing all needed variables.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the pooled ROC curve.}
  \item{B}{An integer value specifying the number of Bayesian bootstrap resamples. By default 5000.}
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
The method implemented in this function makes use of the equivalence (see Gu et al., 2008)
\deqn{ROC(p) = 1 - F_{D}\{F_{\bar{D}}^{-1}(1-p)\} = 1 - Pr(1 - F_{\bar{D}}(Y_D) \leq p),}
and estimates both \eqn{F_{\bar{D}}} and the outer probability using the Bayesian bootstrap resampling distribution. 

Regarding the area under the curve, we note that
\deqn{AUC = \int_{0}^{1}ROC(p)dp = 1 - E\{U_D\},}
where \eqn{U_D = 1 - F_{\bar{D}}(Y_D)}. In our implementation, the expectation is computed using the Bayesian bootstrap (using the same weights as those used to estimate the pooled ROC). As far as the partial area under the curve is concerned, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is 
\deqn{pAUC_{FPF}(u_1)=\int_0^{u_1} ROC(p)dp = u_1 - E\{U_{D,u_1}\},} 
where \eqn{U_{D,u_1} = \min\{u_1, 1 - F_{\bar{D}}(Y_D)\}}. Again, the expectation is computed using the Bayesian bootstrap. The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as 
\deqn{pAUC_{TPF}(u_2)=\int_{u_2}^{1}ROC_{TNF}(p)dp,} 
where \eqn{ROC_{TNF}(p)} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as \eqn{ROC_{TNF}(p) = F_{\bar{D}}\{F_{D}^{-1}(1-p)\}}. Thus
\deqn{pAUC_{TPF}(u_2)=\int_{u_2}^{1}ROC_{TNF}(p)dp = E\{U_{\bar{D}, u_2} - u_2\},} 
where \eqn{U_{\bar{D}, u_2} = \max\{u_2, 1 - F_{D}(Y_{\bar{D}})\}}, and the expectation is computed using the Bayesian bootstrap. The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{marker}{A list with the diagnostic test outcomes in the healthy (h) and diseased (d) groups.}
\item{missing.ind}{A logical value indicating whether missing values occur.}
\item{p}{Set of false positive fractions (FPF) at which the pooled ROC curve has been estimated.}
\item{ROC}{Estimated pooled ROC curve, and corresponding 95\% pointwise credible band.}
\item{AUC}{Estimated pooled AUC, and corresponding 95\% credible interval.}
\item{pAUC}{If computed, estimated partial area under the pooled ROC curve (posterior mean) and 95\% credible interval. Note that the returned values are normalised, so that the maximum value is one (see more on Details).}
\item{weights}{list with the Dirichlet weights (involved in the estimation) in the healthy (h) and diseased (d) groups. These are matrices of dimension n0 x B and n1 x B, where n0 is the number of healthy individuals and n1 is the number of diseased individuals. }
}
\references{
Gu, J., Ghosal, S., and Roy, A. (2008). Bayesian bootstrap estimation of ROC curve. Statistics in Medicine, \bold{27}, 5407--5420.
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
m0_BB <- pooledROC.BB(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 5000)

summary(m0_BB)

plot(m0_BB)
}
\dontshow{
m0_BB <- pooledROC.BB(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 500)

summary(m0_BB)

plot(m0_BB)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
