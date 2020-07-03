\name{cROC.kernel}
\alias{cROC.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric kernel-based estimation of the covariate-specific ROC curve (cROC).
}
\description{
This function estimates the covariate-specific ROC curve (cROC) using the nonparametric kernel-based method proposed by Rodriguez-Alvarez et al. (2011). The method, as it stands now, can only deal with one continuous covariate. 
}
\usage{
cROC.kernel(marker, covariate, group, tag.h, 
  bw = c("LS", "AIC"), regtype = c("LC", "LL"), 
  data, newdata, pauc = pauccontrol(),  
  p = seq(0, 1, l = 101), B = 1000,
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{A character string with the name of the diagnostic test variable.}
  \item{covariate}{A character string with the name of the continuous covariate.}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{bw}{A character string specifying which method to use to select the bandwidths. AIC specifies expected Kullback-Leibler cross-validation, and LS specifies least-squares cross-validation. Defaults to LS. For details see \code{R}-package \code{\link[np:np-package]{np}}.}
  \item{regtype}{A character string specifying which type of kernel estimator to use for the regression function (see Details). LC specifies a local-constant estimator (Nadaraya-Watson) and LL specifies a local-linear estimator. Defaults to LC. For details see \code{R}-package \code{\link[np:np-package]{np}}.}
  \item{data}{Data frame representing the data and containing all needed variables.}
  \item{newdata}{Optional data frame containing the values of the covariates at which the covariate-specific ROC curve (AUC and pAUC, if computed) will be computed. If not supplied, the function \code{\link{cROCData}} is used to build a default dataset.}
   \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve should be computed, and in case it is computed, , whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-specific ROC curve. This set is also used to compute the area under the covariate-specific ROC curve using Simpson's rule. Thus, the length of the set should be an odd number, and it should be rich enough for an accurate estimation.}
  \item{B}{An integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. The default is 1000.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates the covariate-specific ROC curve (cROC) defined as
\deqn{ROC(p|x) = 1 - F_{D}\{F_{\bar{D}}^{-1}(1-p|x)|x\},}
where
\deqn{F_{D}(y|x) = Pr(Y_{D} \leq y | X_{D} = x ),}
\deqn{F_{\bar{D}}(y|x) = Pr(Y_{\bar{D}} \leq y | X_{\bar{D}} = x).}
Note that, for the sake of clarity, we assume that the covariate of interest is the same in both healthy and diseased populations. In particular, the method implemented in this function estimates \eqn{F_{D}(\cdot|x)} and \eqn{F_{\bar{D}}(\cdot|x)} assuming a nonparametric location-scale regression model for \eqn{Y} in each population separately, i.e., 
\deqn{Y_{D} = \mu_{D}(X_{D}) + \sigma_{D}(X_{D})\varepsilon_{D},}
\deqn{Y_{\bar{D}} = \mu_{\bar{D}}(X_{\bar{D}}) + \sigma_{\bar{D}}(X_{\bar{D}})\varepsilon_{\bar{D}},}
where \eqn{\mu_{D}(x) = E(Y_D | X_D = x)},  \eqn{\mu_{\bar{D}}(x) = E(Y_{\bar{D}} | X_{\bar{D}} = x)} (regression function), \eqn{\sigma^2_{D}(x) = Var(Y_D | X_D = x)}, \eqn{\sigma^2_{\bar{D}}(x) = Var(Y_{\bar{D}} | X_{\bar{D}} = x)} (variance functions), and \eqn{\varepsilon_{D}} and \eqn{\varepsilon_{\bar{D}}} have zero mean, variance one, and distribution functions \eqn{G_{D}} and \eqn{G_{\bar{D}}}, respectively. In this case, the covariate-specific ROC curve can be expressed as
\deqn{ROC(p|x) = 1 - G_{D}\{a(\mathbf{x}) + b(\mathbf{x})G_{\bar{D}}^{-1}(1-p)\},}
where \eqn{a(x) = \frac{\mu_{\bar{D}}(x) - \mu_{D}(x)}{\sigma_{D}(x)}}, \eqn{b(x) = \frac{\sigma_{\bar{D}}(x)}{\sigma_{D}(x)}}, and \eqn{G_{D}} and \eqn{G_{\bar{D}}} are the distribution functions of \eqn{\varepsilon_{D}} and \eqn{\varepsilon_{\bar{D}}}, respectively.
By default, for both the healthy and diseased population, both the regression and variance functions are estimated using the Nadaraya-Watson estimator (LC), and the bandwidth are selected using least-squares cross-validation (LS). Implementation relies on the \code{R}-package \code{\link[np:np-package]{np}}. No assumptions are made about \eqn{G_{D}} and \eqn{G_{\bar{D}}}, which are empirically estimated on the basis of standardised residuals. 

The covariate-specific area under the curve is
\deqn{AUC(\mathbf{x})=\int_{0}^{1}ROC(p|\mathbf{x})dp,}
and is computed numerically (using Simpson's rule). With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is
\deqn{pAUC_{FPF}(u_1|\mathbf{x})=\int_0^{u_1} ROC(p|\mathbf{x})dp,}
where again the integral is approximated numerically (Simpson's rule). The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1|\mathbf{x})/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as
\deqn{pAUC_{TPF}(u_2|\mathbf{x})=\int_{u_2}^{1}ROC_{TNF}(p|\mathbf{x})dp,}
where \eqn{ROC_{TNF}(p|\mathbf{x})} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as \eqn{ROC_{TNF}(p|\mathbf{x}) = F_{\bar{D}}\{F_{D}^{-1}(1-p|\mathbf{x})|\mathbf{x}\}=G_{\bar{D}}\{\frac{\mu_{D}(x)-\mu_{\bar{D}}(x)}{\sigma_{\bar{D}}(x)}+G_{D}^{-1}(1-p)\frac{\sigma_{D}(x)}{\sigma_{\bar{D}}(x)}\}.} Again, the computation of the integral is done via Simpson's rule. The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2|\mathbf{x})/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{newdata}{A data frame containing the values of the covariates at which the covariate-specific ROC curve (AUC and pAUC, if required) was computed.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-specific ROC curve has been estimated.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{covariate}{The value of the argument \code{covariate} used in the call.}
\item{ROC}{Estimated covariate-specific ROC curve (AROC), and 95\% pointwise confidence band (if computed).}
\item{AUC}{Estimated area under the covariate-specific ROC curve, and 95\% confidence interval (if computed).}
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve and 95\%  confidence interval (if computed). Note that the returned values are normalised, so that the maximum value is one.}
\item{fit}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component of the list contains the following information: (1) \code{bw.mean}: An object of class \code{\link[np:np.regression.bw]{npregbw}} with the selected bandwidth for the nonparametric regression function. For further details, see \code{R}-package \code{\link[np:np-package]{np}}. (2) \code{bw.var}: An object of class \code{\link[np:np.regression.bw]{npregbw}} with the selected bandwidth for the nonparametric variance function. For further details, see  \code{R}-package \code{\link[np:np-package]{np}}. (3) \code{fit.mean}: An object of class \code{\link[np:np.regression]{npreg}} with the nonparametric regression function estimate. For further details, see \code{R}-package \code{\link[np:np-package]{np}}. (4) \code{fit.var}: An object of class \code{\link[np:np.regression]{npreg}} with the nonparametric variance function estimate. For further details, see \code{R}-package \code{\link[np:np-package]{np}}.}
}
\references{
Hayfield, T., and Racine, J. S.(2008). Nonparametric Econometrics: The np Package. Journal of Statistical Software 27(5). URL http://www.jstatsoft.org/v27/i05/.

Rodriguez-Alvarez, M. X., Roca-Pardinas, J., and Cadarso-Suarez, C. (2011). ROC curve and covariates: extending induced methodology to the non-parametric framework. Statistics and Computing, \bold{21}, 483--499.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{AROC.bnp}}, \code{\link{AROC.sp}}, \code{\link{AROC.kernel}}, \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}}, \code{\link{pooledROC.dpm}}, \code{\link{cROC.kernel}} or \code{\link{cROC.sp}}.
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
cROC_kernel <- cROC.kernel(marker = "l_marker1",
               covariate = "age",
               group = "status", 
               tag.h = 0,
               data = newpsa, 
               bw = "LS",
               regtype = "LC",
               p = seq(0, 1, len = 101),
               pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
               B = 500)

plot(cROC_kernel)

summary(cROC_kernel )
}
\dontshow{
cROC_kernel <- cROC.kernel(marker = "l_marker1",
               covariate = "age",
               group = "status", 
               tag.h = 0,
               data = newpsa, 
               bw = "LS",
               regtype = "LC",
               p = seq(0, 1, len = 101),
               pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
               B = 0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
