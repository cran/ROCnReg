\name{AROC.kernel}
\alias{AROC.kernel}
\title{
Nonparametric kernel-based estimation of the covariate-adjusted ROC curve (AROC).
}
\description{
This function estimates the covariate-adjusted ROC curve (AROC) using the nonparametric kernel-based method proposed by Rodriguez-Alvarez et al. (2011). The method, as it stands now, can only deal with one continuous covariate. 
}
\usage{
AROC.kernel(marker, covariate, group, tag.h, 
    bw = c("LS", "AIC"), 
    regtype = c("LC", "LL"),
    pauc = pauccontrol(), 
    data, p = seq(0, 1, l = 101), B = 1000,
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
   \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve (pAAUC) should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{data}{A data frame representing the data and containing all needed variables.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-adjusted ROC curve. This set is also used to compute the area under the covariate-adjusted ROC curve (AAUC) using Simpson's rule. Thus, the length of the set should be an odd number and it should be rich enough for an accurate estimation.}
  \item{B}{An integer value specifying the number of bootstrap resamples for the construction of the confidence intervals. The default is 1000.}
 \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
}
\details{
Estimates the covariate-adjusted ROC curve (AROC) defined as
\deqn{AROC\left(p\right) = Pr\{1 - F_{\bar{D}}(Y_D | X_{D}) \leq p\},}
where \eqn{F_{\bar{D}}(y|x) = Pr\{Y_{\bar{D}} \leq y | X_{\bar{D}} = x\}}. In particular, the method implemented in this function estimates the outer probability empirically (see Janes and Pepe, 2009) and \eqn{F_{\bar{D}}(y|x)} is estimated assuming a nonparametric location-scale regression model for \eqn{Y_{\bar{D}}}, i.e., 
\deqn{Y_{\bar{D}} = \mu_{\bar{D}}(X_{\bar{D}}) + \sigma_{\bar{D}}(X_{\bar{D}})\varepsilon_{\bar{D}},}
where \eqn{\mu_{\bar{D}}(x) = E(Y_{\bar{D}} | X_{\bar{D}} = x)} is the regression funcion, \eqn{\sigma^2_{\bar{D}}(x) = Var(Y_{\bar{D}} | X_{\bar{D}} = x)} is the variance function, and \eqn{\varepsilon_{\bar{D}}} has zero mean, variance one, and distribution function \eqn{G_{\bar{D}}}. As a consequence,
\deqn{F_{\bar{D}}(y | x) = G_{\bar{D}}\left(\frac{y - \mu_{\bar{D}}(x)}{\sigma_{\bar{D}}(x)}\right).} 
By default, both the regression and variance functions are estimated using the Nadaraya-Watson estimator (LC), and the bandwidths are selected using least-squares cross-validation (LS). Implementation relies on the \code{R}-package \code{\link[np:np-package]{np}}. No assumption is made about \eqn{G_{\bar{D}}}, which is empirically estimated on the basis of the standardised residuals. 

The area under the AROC curve is
\deqn{AAUC=\int_0^1 AROC(p)dp,}
and there exists a closed-form estimator. With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is
\deqn{pAAUC_{FPF}(u_1)=\int_0^{u_1} AROC(p)dp,}
where again there exists a closed-form estimator. The returned value is the normalised pAAUC, \eqn{pAAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as 
\deqn{pAAUC_{TPF}(u_2)=\int_{AROC^{-1}(u_2)}^{1}AROC(p)dp-\{1-AROC^{-1}(u_2)\}\times u_2.}
Here, the computation of the integral is done numerically. The returned value is the normalised pAAUC, \eqn{pAAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-adjusted ROC curve has been estimated.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{covariate}{The value of the argument \code{covariate} used in the call.}
\item{ROC}{Estimated covariate-adjusted ROC curve (AROC), and 95\% pointwise confidence band (if computed).}
\item{AUC}{Estimated area under the covariate-adjusted ROC curve (AAUC), and 95\% confidence interval (if computed).}
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve (pAAUC) and 95\%  confidence interval (if computed). Note that the returned values are normalised, so that the maximum value is one.}
\item{fit}{List with the following components: (1) \code{bw.mean}: An object of class \code{\link[np:np.regression.bw]{npregbw}} with the selected bandwidth for the nonparametric regression function. For further details, see  \code{R}-package \code{\link[np:np-package]{np}}. (2) \code{bw.var}: An object of class \code{\link[np:np.regression.bw]{npregbw}} with the selected bandwidth for the nonparametric variance function. For further details, see  \code{R}-package \code{\link[np:np-package]{np}}. (3) \code{fit.mean}: An object of class \code{\link[np:np.regression]{npreg}} with the nonparametric regression function estimate. For further details, see  \code{R}-package \code{\link[np:np-package]{np}}. (4) \code{fit.var}: An object of class \code{\link[np:np.regression]{npreg}} with the nonparametric variance function estimate. For further details, see  \code{R}-package \code{\link[np:np-package]{np}}.}
}
\references{
Hayfield, T., and Racine, J. S. (2008). Nonparametric Econometrics: The np Package. Journal of Statistical Software 27(5). URL http://www.jstatsoft.org/v27/i05/.

Inacio de Carvalho, V., and Rodriguez-Alvarez, M. X. (2018). Bayesian nonparametric inference for the covariate-adjusted ROC curve. arXiv preprint arXiv:1806.00473.

Janes, H., and Pepe, M.S. (2009). Adjusting for covariate effects on classification accuracy using the covariate-adjusted receiver operating characteristic curve. Biometrika, \bold{96}, 371--382.

Rodriguez-Alvarez, M. X., Roca-Pardinas, J., and Cadarso-Suarez, C. (2011). ROC curve and covariates: extending induced methodology to the non-parametric framework. Statistics and Computing, \bold{21}, 483--499.
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
m2 <- AROC.kernel(marker = "l_marker1", 
covariate = "age",
group = "status", 
tag.h = 0,
data = newpsa, 
bw = "LS",
regtype = "LC",
pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
B = 500)

summary(m2)

plot(m2)
}
\dontshow{
m2 <- AROC.kernel(marker = "l_marker1", 
covariate = "age",
group = "status", 
tag.h = 0,
data = newpsa, 
bw = "LS",
regtype = "LC",
pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
B = 0)

summary(m2)

plot(m2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
