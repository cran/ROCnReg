\name{AROC.bnp}
\alias{AROC.bnp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Bayesian inference of the covariate-adjusted ROC curve (AROC).
}
\description{
This function estimates the covariate-adjusted ROC curve (AROC) using the nonparametric Bayesian approach proposed by Inacio de Carvalho and Rodriguez-Alvarez (2018).
}
\usage{
AROC.bnp(formula.h, group, tag.h, data, standardise = TRUE, 
  p = seq(0, 1, l = 101), compute.lpml = FALSE, compute.WAIC = FALSE, 
  compute.DIC = FALSE, pauc = pauccontrol(), density = densitycontrol.aroc(),
  prior.h = priorcontrol.bnp(), mcmc = mcmccontrol(),
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.h}{A \code{\link[stats]{formula}} object specifying the regression function associated to each component of the single-weights linear dependent Dirichlet process mixture model used to estimate the conditional distribution function of the diagnostic test outcome in the healthy population. Regarding the modelling of continuous covariates, both linear and nonlinear effects are allowed, with nonlinear effects being modelled through B-spline basis expansions (see Note).}
  \item{group}{A character string with the name of the variable that distinguishes healthy/nondiseased from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{ A data frame representing the data and containing all needed variables.}
  \item{standardise}{A logical value. If TRUE both the test outcomes and the continuous covariates assumed to have a linear effect are standardised (i.e., the resulting variables have mean zero and standard deviation of one). The default is TRUE.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-adjusted ROC curve.}
  \item{compute.lpml}{A logical value. If TRUE, the log pseudo marginal likelihood (LPML, Geisser and Eddy, 1979) and the conditional predictive ordinates (CPO) are computed.}
  \item{compute.WAIC}{A logical value. If TRUE, the widely applicable information criterion (WAIC, Gelman et al., 2014; Watanabe, 2010) is computed.}
  \item{compute.DIC}{A logical value. If TRUE, the deviance information criterion (DIC)(Celeux et al., 2006, Spiegelhalter et al., 2002) is computed.}
  \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve (pAAUC) should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{density}{A list of control values to replace the default values returned by the function \code{\link{densitycontrol.aroc}}. This argument is used to indicate whether the conditional densities of the marker in the healthy population should be computed, and in case it is to be computed, at which grid of test outcomes the conditional densities should be evaluated, and at which covariate values they should be predicted.}
  \item{prior.h}{A list of control values to replace the default values returned by the function \code{\link{priorcontrol.bnp}}. See \code{\link{priorcontrol.bnp}} for details.}
  \item{mcmc}{A list of control values to replace the default values returned by the function \code{\link{mcmccontrol}}. See \code{\link{mcmccontrol}} for details.}
  \item{parallel}{A characters string with the type of parallel operation: either "no" (default), "multicore" (not available on Windows) or "snow".}
  \item{ncpus}{An integer with the number of processes to be used in parallel operation. Defaults to 1.}
  \item{cl}{An object inheriting from class \code{cluster} (from the \code{parallel} package), specifying an optional parallel or snow cluster if parallel = "snow". If not supplied, a cluster on the local machine is created for the duration of the call.}
} 
 
\details{
Estimates the covariate-adjusted ROC curve (AROC) defined as
\deqn{AROC\left(p\right) = Pr\{1 - F_{\bar{D}}(Y_D | \mathbf{X}_{D}) \leq p\},}
where \eqn{F_{\bar{D}}(\cdot|\mathbf{X}_{\bar{D}})} denotes the distribution function of \eqn{Y_{\bar{D}}} conditional on the vector of covariates \eqn{\mathbf{X}_{\bar{D}}}. 

The method implemented in this function combines a single-weights linear dependent Dirichlet process mixture model (De Iorio et al., 2009) to estimate \eqn{F_{\bar{D}}(\cdot|\mathbf{X}_{\bar{D}})} and the Bayesian bootstrap (Rubin, 1981) to estimate the outside probability. More precisely, and letting \eqn{\{(\mathbf{x}_{\bar{D}i},y_{\bar{D}i})\}_{i=1}^{n_{\bar{D}}}} be a random sample from the nondiseased population, our postulated model for the conditional distribution function takes the following form 
\deqn{F_{\bar{D}}(y_{\bar{D}i}|\mathbf{X}_{\bar{D}}=\mathbf{x}_{\bar{D}i}) = \sum_{l=1}^{L}\omega_l\Phi(y_{\bar{D}i}\mid\mu_{l}(\mathbf{x}_{\bar{D}i}),\sigma_l^2),}
where \eqn{\Phi(y|\mu, \sigma^2)} denotes the cumulative distribution function of the normal distribution, evaluated at \eqn{y}, with mean \eqn{\mu} and variance \eqn{\sigma^2}. The regression function \eqn{\mu_{l}(\mathbf{x}_{\bar{D}i})} can incorportate both linear and nonlinear (through B-splines) effects of continuous covariates, categorical covariates (factors) as well as interactions. Interactions between categorical and (nonlinear) continuous covariates are also allowed (factor-by curve interactions). For the sake of simplicity we write \eqn{\mu_{l}(\mathbf{x}_{\bar{D}i}) = \mathbf{z}_{\bar{D}i}^{T}\mathbf{\beta}_l} (\eqn{l=1,...,L}), where \eqn{\mathbf{z}_{\bar{D}i}} is the \eqn{i}th column of the design matrix (possibly containing a basis representation of some/all continuous covariates). Here \eqn{L} is a pre-specified upper bound on the number of mixture components. The \eqn{\omega_l}'s result from a truncated version of the stick-breaking construction (\eqn{\omega_1=v_1}; \eqn{\omega_l=v_l\prod_{r<l}(1-v_r)}, \eqn{l=2,\ldots,L}; \eqn{v_1,\ldots,v_{L-1}\sim} Beta \eqn{(1,\alpha)}; \eqn{v_L=1}, \eqn{\alpha \sim \Gamma(a_{\alpha},b_{\alpha})}), \eqn{\mathbf{\beta}_l\sim N_{Q}(\mathbf{m},\mathbf{S})}, and \eqn{\sigma_l^{-2}\sim\Gamma(a,b)}. It is further assumed that \eqn{\mathbf{m} \sim N_{Q}(\mathbf{m}_0,\mathbf{S}_0)} and \eqn{\mathbf{S}^{-1}\sim W(\nu,(\nu\Psi)^{-1})}. Here \eqn{\Gamma(a,b)} denotes a Gamma distribution with shape parameter \eqn{a} and rate parameter \eqn{b}, \eqn{W(\nu,(\nu\Psi)^{-1})} denotes a Wishart distribution with \eqn{\nu} degrees of freedom and expectation \eqn{\Psi^{-1}}, and \eqn{Q} denotes the dimension of the vector \eqn{\mathbf{z}_{\bar{D}i}}. It is worth mentioning that when \eqn{L=1}, the model for the conditional distribution of the test outcomes (in the healthy population) reduces to a normal regression model (where continuous covariates effects are modelled either parametrically or nonparametrically). For a detailed description, we refer to Inacio de Carvalho and Rodriguez-Alvarez (2018).

Regarding the area under the curve, we note that
\deqn{AAUC = \int_{0}^{1}AROC(p)dp = 1 - E\{U_D\},}
where \eqn{U_D = 1 - F_{\bar{D}}(Y_D |\mathbf{X}_D)}. In our implementation, the expectation is computed using the Bayesian bootstrap (using the same weights as those used to estimate the AROC, see Inacio de Carvalho and Rodriguez-Alvarez (2018) for details). As far as the partial area under the curve is concerned, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is 
\deqn{pAAUC_{FPF}(u_1)=\int_0^{u_1} AROC(p)dp = u_1 - E\{U_{D,u_1}\},} 
where \eqn{U_{D,u_1} = min\{u_1, 1 - F_{\bar{D}}(Y_D |\mathbf{X}_D)\}}. Again, the expectation is computed using the Bayesian bootstrap. The returned value is the normalised pAAUC, \eqn{pAAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as 
\deqn{pAAUC_{TPF}(u_2)=\int_{AROC^{-1}(u_2)}^{1}AROC(p)dp-\{1-AROC^{-1}(u_2)\}\times u_2.}
Here, the computation of the integral is done numerically. The returned value is the normalised pAAUC, \eqn{pAAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

Finally, it is important referring that with respect to the computation of the DIC, when \eqn{L=1}, it is computed as in Spiegelhalter et al. (2002), and when \eqn{L>1}, DIC3 as described in Celeux et al. (2006) is computed. Also, for the computation of the conditional predictive ordinates (CPO) we follow the stable version proposed by Gelman et al. (2014).

}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-adjusted ROC curve (AROC) has been estimated.}
\item{prior}{A list returning the hyperparameter values.}
\item{ROC}{Estimated covariate-adjusted ROC curve (AROC) (posterior mean) and 95\% pointwise credible band.}
\item{AUC}{Estimated area under the covariate-adjusted ROC curve (AAUC) (posterior mean), and 95\% credible interval.}
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve (pAAUC) (posterior mean) and 95\% credible
interval. Note that the returned values are normalised, so that the maximum value is one (see more on Details).}
\item{newdata}{If \code{compute} is set to \code{TRUE} in the argument \code{density}, a data frame containing the values of the covariates at which the regression function and conditional densities were computed (see below).}
\item{reg.fun.h}{If \code{compute} is set to \code{TRUE} in the argument \code{density}, a data frame containing the predicted regression function (posterior mean) and 95\% pointwise credible band.}
\item{dens}{If \code{compute} is set to \code{TRUE} in the argument \code{density}, a list with two components (only for the healthy population): \code{grid} (grid of test outcomes where the densities were evaluated) and \code{dens} (MCMC realisations of the corresponding conditional densities).}
\item{lpml}{If computed, a list with two components: the log pseudo marginal likelihood (LPML) and the conditional predictive ordinates (CPO).}
\item{WAIC}{If computed, widely applicable information criterion (WAIC) and associated complexity penalty (pW).}
\item{DIC}{If computed, deviance information criterion (DIC) and associated complexity penalty (pD).}
\item{fit}{Results of the fitting process. A list with the following components: (1) \code{formula}: the value of the argument \code{formula.h} used in the call. (2) \code{mm}: information needed to construct the model design matrix associated with the single weights linear dependent Dirichlet process mixture model. (3) \code{beta}: array of dimension \code{nsave}x\code{L}x\code{Q} with the sampled regression coefficients. (4) \code{sd}: matrix of dimension \code{nsave}x\code{L} with the sampled variances. (4) \code{probs}: matrix of dimension \code{nsave}x\code{L} with the sampled components' weights. Here, \code{nsave} is the number of Gibbs sampler iterations saved, \code{L} is the upper bound on the number of mixture components, and \code{Q} is the dimension of vector \eqn{\mathbf{z}_{\bar{D}}} (see also Details).}
\item{data_model}{List with the data used in the fit: observed diagnostic test outcome and design matrices, separately for the healthy and diseased groups.}
}
\references{
Celeux, G., Forbes, F., Robert C. P., and Titerrington, D. M. (2006). Deviance information criteria for missing data models. Bayesian Analysis, \bold{1}, 651--674.

De Iorio, M., Johnson, W. O., Muller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. Biometrics, \bold{65}, 762--775.

Geisser, S. and Eddy, W.F. (1979) A Predictive Approach to Model Selection, Journal of the American Statistical Association, \bold{74}, 153--160.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and Rubin, D.B. (2014). Bayesian Data Analysis, 3rd ed. CRC Press: Boca Raton, FL.

Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Statistics and Computing, \bold{24}, 997--1010.

Inacio de Carvalho, V., and Rodriguez-Alvarez, M. X. (2018). Bayesian nonparametric inference for the covariate-adjusted ROC curve. arXiv preprint arXiv:1806.00473.

Rubin, D. B. (1981). The Bayesian bootstrap. The Annals of Statistics, \bold{9}, 130--134.

Speigelhalter, D. J., Best, N. G., Carlin, B. P., and van der Linde, A. (2002). Bayesian measures of model comparison and fit. Journal of the Royal Statistical Society, Ser. B, \bold{64}, 583--639.

Watanabe, S. (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular Learning Theory. Journal of Machine Learning Research, \bold{11}, 3571--3594.
}
%\author{
%%  ~~who you are~~
%}
\note{
The input argument \code{formula.h} is similar to that used for the \code{\link{glm}} function, except that flexible specifications can be added by means of the function \code{f()}. For instance, specification \eqn{y \sim x1 + f(x2, K = 3)} would assume a linear effect of \code{x1} (if \code{x1} continuous) and the effect of \code{x2} would be modeled using B-splines basis functions. The argument \code{K = 3} indicates that \code{3} internal knots will be used, with the quantiles of \code{x2} used for their location. Categorical variables (factors) can be also incorporated, as well as interaction terms. For example, to include the factor-by-curve interaction between \code{age} and \code{gender} we need to specify, e.g., \eqn{y \sim gender + f(age, by = gender, K = c(3, 5))}. Note that, in this case, the number of knots can be different for each level of the factor. The order of the vector \code{K} of knots should match the levels of the factor.}

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
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              compute.lpml = TRUE,
              compute.WAIC = TRUE,
              compute.DIC = TRUE,
              pauc = pauccontrol(compute = TRUE, focus = "FPF", value = 0.5),
              density = densitycontrol.aroc(compute = TRUE, grid.h = NA, newdata = NA),
              prior.h = priorcontrol.bnp(m0 = rep(0, 4), S0 = 10*diag(4), nu = 6, Psi = diag(4),
              a = 2, b = 2, aalpha = 2, balpha = 2, L =10),
              mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1))

summary(AROC_bnp)

plot(AROC_bnp)
}
\dontshow{

AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
              group = "status", 
              tag.h = 0, 
              data = newpsa,
              standardise = TRUE,
              p = seq(0,1,l=101),
              compute.lpml = FALSE,
              compute.WAIC = FALSE,
              compute.DIC = FALSE,
              pauc = pauccontrol(compute = FALSE, focus = "FPF", value = 0.5),
              density = densitycontrol.aroc(compute = FALSE, grid.h = NA, newdata = NA),
              prior.h = priorcontrol.bnp(m0 = rep(0, 4), S0 = 10*diag(4), nu = 6, Psi = diag(4),
              a = 2, b = 2, aalpha = 2, balpha = 2, L =10),
              mcmc = mcmccontrol(nsave = 50, nburn = 10, nskip = 1))

summary(AROC_bnp)

plot(AROC_bnp)

}
}