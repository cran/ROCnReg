\name{cROC.sp}
\alias{cROC.sp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric and semiparametric frequentist inference of the covariate-specific ROC curve (cROC).
}
\description{
This function estimates the covariate-specific ROC curve (cROC) using the parametric approach proposed by Faraggi (2003) and the semiparametric approach proposed by Pepe (1998).
}
\usage{
cROC.sp(formula.h, formula.d, group, tag.h, data, 
  newdata, est.cdf = c("normal", "empirical"),
  pauc = pauccontrol(), p = seq(0, 1, l = 101), B = 1000,
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.h}{A \code{\link[stats]{formula}} object specifying the location regression model to be fitted in the healthy population (see Details).}
  \item{formula.d}{A \code{\link[stats]{formula}} object specifying the location regression model to be fitted in the diseased population (see Details).}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying the healthy individuals in the variable \code{group}.}
  \item{data}{Data frame representing the data and containing all needed variables.}
  \item{newdata}{Optional data frame containing the values of the covariates at which the covariate-specific ROC curve (AUC and pAUC, if required) will be computed. If not supplied, the function \code{\link{cROCData}} is used to build a default dataset.}
  \item{est.cdf}{A character string. It indicates how the conditional distribution functions of the diagnostic test in healthy and diseased populations are estimated. Options are "normal" and "empirical" (see Details). The default is "normal".}
  \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-specific ROC curve. This set is also used to compute the area under the covariate-specific ROC curve using Simpson's rule. Thus, the length of the set should be an odd number, and it should be rich enough for an accurate estimation.}
  \item{B}{An integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. By default 1000.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates the covariate-specific ROC curve (cROC) defined as
\deqn{ROC(p|\mathbf{x}) = 1 - F_{D}\{F_{\bar{D}}^{-1}(1-p|\mathbf{x})|\mathbf{x}\},}
where
\deqn{F_{D}(y|\mathbf{x}) = Pr(Y_{D} \leq y | \mathbf{X}_{D} = \mathbf{x}),}
\deqn{F_{\bar{D}}(y|\mathbf{x}) = Pr(Y_{\bar{D}} \leq y | \mathbf{X}_{\bar{D}} = \mathbf{x}).}
Note that, for the sake of clarity, we assume that the covariates of interest are the same in both healthy and diseased populations. In particular, the method implemented in this function estimates \eqn{F_{D}(\cdot|\mathbf{x})} and \eqn{F_{\bar{D}}(\cdot|\mathbf{x})} assuming a (semiparametric) location regression model for \eqn{Y} in each population separately, i.e., 
\deqn{Y_{D} = \mathbf{X}_{D}^{T}\mathbf{\beta}_{D} + \sigma_{D}\varepsilon_{D},}
\deqn{Y_{\bar{D}} = \mathbf{X}_{\bar{D}}^{T}\mathbf{\beta}_{\bar{D}} + \sigma_{\bar{D}}\varepsilon_{\bar{D}},}
such that the covariate-specific ROC curve can be expressed as
\deqn{ROC(p|\mathbf{x}) = 1 - G_{D}\{a(\mathbf{x}) + b G_{\bar{D}}^{-1}(1-p)\},}
where \eqn{a(\mathbf{x}) = \mathbf{x}^{T}\frac{\mathbf{\beta}_{\bar{D}} - \mathbf{\beta}_{D}}{\sigma_{D}}}, \eqn{b = \frac{\sigma_{\bar{D}}}{\sigma_{D}}}, and \eqn{G_{D}} and \eqn{G_{\bar{D}}} are the distribution functions of \eqn{\varepsilon_{D}} and \eqn{\varepsilon_{\bar{D}}}, respectively. In line with the assumptions made about the distributions of \eqn{\varepsilon_{D}} and \eqn{\varepsilon_{\bar{D}}}, estimators will be referred to as: (a) "normal", where Gaussian errors are assumed, i.e., \eqn{G_{D}(y) = G_{\bar{D}}(y) = \Phi(y)} (Faraggi, 2003); and, (b) "empirical", where no assumptios are made about the distribution (in this case, \eqn{G_{D}} and \eqn{G_{\bar{D}}} are empirically estimated on the basis of standardised residuals (Pepe, 1998)).

The covariate-specific area under the curve is
\deqn{AUC(\mathbf{x})=\int_{0}^{1}ROC(p|\mathbf{x})dp.}
When Gaussian errors are assumed, there is a closed-form expression for the covariate-specific AUC, which is used in the package. In contrast, when no assumptios are made about the distributionis of the errors, the integral is computed numerically using Simpson's rule. With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is
\deqn{pAUC_{FPF}(u_1|\mathbf{x})=\int_0^{u_1} ROC(p|\mathbf{x})dp.}
Again, when Gaussian errors are assumed, there is a closed-form expression (Hillis and Metz, 2012). Otherwise, the integral is approximated numerically (Simpson's rule). The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1|\mathbf{x})/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as
\deqn{pAUC_{TPF}(u_2|\mathbf{x})=\int_{u_2}^{1}ROC_{TNF}(p|\mathbf{x})dp,}
where \eqn{ROC_{TNF}(p|\mathbf{x})} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as\eqn{ROC_{TNF}(p|\mathbf{x}) = F_{\bar{D}}\{F_{D}^{-1}(1-p|\mathbf{x})|\mathbf{x}\}=G_{\bar{D}}\{\frac{\mu_{D}(\mathbf{x})-\mu_{\bar{D}}(\mathbf{x})}{\sigma_{\bar{D}}(\mathbf{x})}+G_{D}^{-1}(1-p)\frac{\sigma_{D}(\mathbf{x})}{\sigma_{\bar{D}}(\mathbf{x})}\}.} Again, when Gaussian errors are assumed, there is a closed-form expression (Hillis and Metz, 2012). Otherwise, the integral is approximated numerically (Simpson's rule). The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2|\mathbf{x})/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

}
\value{As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{newdata}{A data frame containing the values of the covariates at which the covariate-specific ROC curve (AUC and pAUC, if required) was computed.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{formula}{Named list of length two with the value of the arguments \code{formula.h} and \code{formula.d} used in the call.}
\item{est.cdf}{The value of the argument \code{est.cdf} used in the call.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-specific ROC curves have been estimated.}
\item{ROC}{Estimated covariate-specific ROC curve, and 95\% pointwise confidence intervals (if computed).}
\item{AUC}{Estimated area under the covariate-specific ROC curve, and 95\% confidence interval (if computed). }
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve and 95\%  confidence interval (if computed). Note that the returned values are normalised, so that the maximum value is one.}
\item{fit}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component contains an object of class \code{\link[stats]{lm}} with the fitted regression model.}
\item{coeff}{Estimated regression coefficients (and 95\% confidence interval if \code{B} greater than zero) from the fit of the linear model in the healthy and diseased population, as specified in \code{formula.h} and \code{formula.d}, respectively.}
}
\references{
Faraggi, D. (2003). Adjusting receiver operating characteristic curves and related indices for covariates. \emph{The Statistician} \bold{52}, 179--192.

Hillis, S. L. and Metz, C.E. (2012). An Analytic Expression for the Binormal Partial Area under the ROC Curve. Academic Radiology, \bold{19}, 1491--1498.

Pepe, M.S. (1998). Three approaches to regression analysis of receiver operating characteristic curves for continuous test results. \emph{Biometrics} \bold{54}, 124--135.

}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

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

# Covariate for prediction
agep <- seq(min(newpsa$age), max(newpsa$age), length = 50)
df.pred <- data.frame(age = agep)

\donttest{
cROC_sp_normal <- cROC.sp(formula.h = l_marker1 ~ age,
                          formula.d = l_marker1 ~ age,
                          group = "status", 
                          tag.h = 0,
                          data = newpsa,
                          newdata = df.pred,
                          est.cdf = "normal",
                          pauc = list(compute = TRUE, value = 0.5, focus = "FPF"),
                          p = seq(0, 1, l = 101), 
                          B = 500)
summary(cROC_sp_normal)

plot(cROC_sp_normal)
}
\dontshow{
cROC_sp_normal <- cROC.sp(formula.h = l_marker1 ~ age,
                          formula.d = l_marker1 ~ age,
                          group = "status", 
                          tag.h = 0,
                          data = newpsa,
                          newdata = df.pred,
                          est.cdf = "normal",
                          pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
                          p = seq(0, 1, l = 101), 
                          B = 10)
summary(cROC_sp_normal)

plot(cROC_sp_normal)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
