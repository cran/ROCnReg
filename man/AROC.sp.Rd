\name{AROC.sp}
\alias{AROC.sp}
\title{
Semiparametric frequentist inference for the covariate-adjusted ROC curve (AROC).
}
\description{
This function estimates the covariate-adjusted ROC curve (AROC) using the semiparametric approach proposed by Janes and Pepe (2009).
}
\usage{
AROC.sp(formula.h, group, tag.h, data, 
    est.cdf.h = c("normal", "empirical"), pauc = pauccontrol(),
    p = seq(0, 1, l = 101), B = 1000,
  	parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.h}{A \code{\link[stats]{formula}} object specifying the location regression model to be fitted in the healthy population (see Details).}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{A data frame representing the data and containing all needed variables.}
  \item{est.cdf.h}{A character string. It indicates how the conditional distribution function of the diagnostic test in the healthy population is estimated. Options are \code{"normal"} and \code{"empirical"} (see Details). The default is \code{"normal"}.}
 \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve (pAAUC) should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-adjusted ROC curve. This set is also used to compute the area under the covariate-adjusted ROC curve (AAUC) using Simpson's rule. Thus, the length of the set should be an odd number, and it should be rich enough for an accurate estimation.}
  \item{B}{An integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. The default is 1000.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates the covariate-adjusted ROC curve (AROC) defined as
\deqn{AROC\left(p\right) = Pr\{1 - F_{\bar{D}}(Y_D | \mathbf{X}_{D}) \leq p\},}
\deqn{F_{\bar{D}}(y|\mathbf{x}) = Pr\{Y_{\bar{D}} \leq y | \mathbf{X}_{\bar{D}} = \mathbf{x}\}.} The method implemented in this function estimates the outer probability empirically (see Janes and Pepe, 2009) and \eqn{F_{\bar{D}}(\cdot|\mathbf{x})} is estimated assuming a semiparametric location regression model for \eqn{Y_{\bar{D}}}, i.e., 
\deqn{Y_{\bar{D}} = \mathbf{X}_{\bar{D}}^{T}\mathbf{\beta}_{\bar{D}} + \sigma_{\bar{D}}\varepsilon_{\bar{D}},}
where \eqn{\varepsilon_{\bar{D}}} has zero mean, variance one, and distribution function \eqn{G_{\bar{D}}}. As a consequence, we have
\deqn{F_{\bar{D}}(y | \mathbf{x}) = G_{\bar{D}}\left(\frac{y-\mathbf{x}^{T}\mathbf{\beta}_{\bar{D}}}{\sigma_{\bar{D}}}\right).}
In line with the assumptions made about the distribution of \eqn{\varepsilon_{\bar{D}}}, estimators will be referred to as: (a) "normal", where a standard Gaussian error is assumed, i.e., \eqn{G_{\bar{D}}(y) = \Phi(y)}; and, (b) "empirical", where no assumption is made about the distribution (in this case, \eqn{G_{\bar{D}}} is empirically estimated on the basis of standardised residuals).

The area under the AROC curve is
\deqn{AAUC=\int_0^1 AROC(p)dp,}
and there exists a closed-form estimator. With regard to the partial area under the AROC curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is
\deqn{pAAUC_{FPF}(u_1)=\int_0^{u_1} AROC(p)dp,}
where again there exists a closed-form estimator. The returned value is the normalised pAAUC, \eqn{pAAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as \deqn{pAAUC_{TPF}(u_2)=\int_{AROC^{-1}(u_2)}^{1}AROC(p)dp-\{1-AROC^{-1}(u_2)\}\times u_2.} Here, the computation of the integral is done numerically. The returned value is the normalised pAAUC, \eqn{pAAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{formula}{The value of the argument \code{formula.h} used in the call.}
\item{est.cdf.h}{The value of the argument \code{est.cdf.h} used in the call.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-adjusted ROC (AROC) curve has been estimated}
\item{ROC}{Estimated covariate-adjusted ROC curve (AROC), and 95\% pointwise confidence bands (if computed)}
\item{AUC}{Estimated area under the covariate-adjusted ROC curve (AAUC), and 95\% confidence intervals (if required).}
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve (pAAUC) and 95\%  confidence interval (if computed). Note that the returned values are normalised, so that the maximum value is one.}
\item{fit}{Object of class \code{\link[stats]{lm}} with the fitted regression model in the healthy population.}
\item{coeff}{Estimated regression coefficients (and 95\% confidence interval if \code{B} greater than zero) from the fit of the linear model in the healthy population, as specified in \code{formula.h}.}
}
\references{
Janes, H., and Pepe, M.S. (2009). Adjusting for covariate effects on classification accuracy using the covariate-adjusted receiver operating characteristic curve. Biometrika, 96(2), 371 - 382.
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
m3 <- AROC.sp(formula.h = l_marker1 ~ age,
group = "status", 
tag.h = 0,
data = newpsa,
est.cdf.h = "normal",
pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
p = seq(0,1,l=101), 
B = 500)

summary(m3)

plot(m3)
}
\dontshow{
m3 <- AROC.sp(formula.h = l_marker1 ~ age,
group = "status", 
tag.h = 0,
data = newpsa,
est.cdf.h = "normal",
pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
p = seq(0,1,l=101), 
B = 10)

summary(m3)

plot(m3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line