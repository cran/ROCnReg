\name{cROC.bnp}
\alias{cROC.bnp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Bayesian inference for the covariate-specific ROC curve (cROC).
}
\description{
This function estimates the covariate-specific ROC curve (cROC) using the nonparametric Bayesian approach proposed by Inacio de Carvalho et al. (2013).
}
\usage{
cROC.bnp(formula.h, formula.d, group, tag.h, data, 
    newdata, standardise = TRUE, p = seq(0, 1, l = 101), 
    compute.lpml = FALSE, compute.WAIC = FALSE, compute.DIC = FALSE, 
    pauc = pauccontrol(), density = densitycontrol(),
    prior.h = priorcontrol.bnp(), prior.d = priorcontrol.bnp(), 
    mcmc = mcmccontrol(),
    parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) 
}
\arguments{
   \item{formula.h}{A \code{\link[stats]{formula}} object specifying the regression function associated to each component of the single-weights linear dependent Dirichlet process mixture of normals model used to estimate the conditional distribution function of the diagnostic test outcome in the healthy population. Regarding the modelling of continuous covariates, both linear and nonlinear effects are allowed, with nonlinear effects being modelled through B-spline basis expansions (see Note).}
  \item{formula.d}{A \code{\link[stats]{formula}} object specifying the regression function associated to each component of the single weights linear dependent Dirichlet process mixture model used to estimate the conditional distribution function of the diagnostic test outcome in the diseased population. Both linear and nonlinear (through the use of smooth functions approximated by linear combinations of cubic B-splines basis functions) covariate effects are allowed (see Note).}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{A data frame representing the data and containing all needed variables.}
  \item{newdata}{Optional data frame containing the values of the covariates at which the covariate-specific ROC curve (AUC and pAUC, if computed) will be computed. If not supplied, the function \code{\link{cROCData}} is used to build a default dataset.}
   \item{standardise}{A logical value. If TRUE both the test outcomes and the continuous covariates assumed to have a linear effect are standardised (i.e., the resulting variables have mean zero and standard deviation of one). The default is TRUE.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the covariate-specific ROC curve.}
  \item{compute.lpml}{A logical value. If TRUE, the log pseudo marginal likelihood (LPML, Geisser and Eddy, 1979) and the conditional predictive ordinates (CPO) are computed.}
  \item{compute.WAIC}{A logical value. If TRUE, the widely applicable information criterion (WAIC, Gelman et al., 2014; Watanabe, 2010) is computed.}
  \item{compute.DIC}{A logical value. If TRUE, the deviance information criterion (DIC)(Celeux et al., 2006, Spiegelhalter et al., 2002) is computed.}
  \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the covariate-adjusted ROC curve should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
  \item{density}{A list of control values to replace the default values returned by the function \code{\link{densitycontrol}}. This argument is used to indicate whether the conditional densities of the marker in the healthy and diseased population should be computed, and in case it is to be computed, at which grid of test outcomes in each of the populations the conditional densities should be evaluated.}
  \item{prior.h}{Hyparameter specification for the healthy population. A list of control values to replace the default values returned by the function \code{\link{priorcontrol.bnp}}. See \code{\link{priorcontrol.bnp}} for details.}
  \item{prior.d}{Hyparameter specification for the diseased population.  A list of control values to replace the default values returned by the function \code{\link{priorcontrol.bnp}}. See \code{\link{priorcontrol.bnp}} for details.}
 \item{mcmc}{A list of control values to replace the default values returned by the function \code{\link{mcmccontrol}}. See \code{\link{mcmccontrol}} for details.}
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
Note that, for the sake of clarity, we assume that the covariates of interest are the same in both the healthy and diseased populations. The method implemented in this function estimates \eqn{F_{D}(\cdot|\mathbf{x})} and \eqn{F_{\bar{D}}(\cdot|\mathbf{x})} by means of a single-weights linear dependent Dirichlet process mixture of normals model (De Iorio et al., 2009). More precisely, and letting \eqn{\{(\mathbf{x}_{\bar{D}i},y_{\bar{D}i})\}_{i=1}^{n_{\bar{D}}}}  and \eqn{\{(\mathbf{x}_{Dj},y_{Dj})\}_{j=1}^{n_{D}}} be two independent random samples from the nondiseased and diseased populations, respectively, our postulated model for the conditional distribution in each group function takes the following form 
\deqn{F_{\bar{D}}(y_{\bar{D}i}|\mathbf{X}_{\bar{D}}=\mathbf{x}_{\bar{D}i}) = \sum_{l=1}^{L_{\bar{D}}}\omega_{l\bar{D}}\Phi(y_{\bar{D}i}\mid\mu_{l\bar{D}}(\mathbf{x}_{\bar{D}i}),\sigma_{l\bar{D}}^2),}
\deqn{F_{D}(y_{Dj}|\mathbf{X}_{D} = \mathbf{x}_{Dj}) = \sum_{l=1}^{L_{D}}\omega_{lD}\Phi(y_{Dj}\mid\mu_{lD}(\mathbf{x}_{\bar{D}i}),\sigma_{lD}^2),}
where \eqn{\Phi(y|\mu, \sigma^2)} denotes the cumulative distribution function of the normal distribution, evaluated at \eqn{y}, with mean \eqn{mu} and variance \eqn{\sigma^2}. The regression function \eqn{\mu_{ld}(\mathbf{x}_{di})} can incorportate both linear and nonlinear (through B-splines) effects of continuous covariates, categorical covariates (factors) as well as interactions. Interactions between categorical and (nonlinear) continuous covariates are also allowed (factor-by curve interactions). For the sake of simplicity we write \eqn{\mu_{ld}(\mathbf{x}_{di}) = \mathbf{z}_{di}^{T}\mathbf{\beta}_{ld}}, where \eqn{\mathbf{z}_{di}} is the \eqn{i}th column of the design matrix (possibly containing a basis representation of some/all continuous covariates), \eqn{d \in \{D, \bar{D}\}}. Here \eqn{L_d} is a pre-specified upper bound on the number of mixture components. The \eqn{\omega_{ld}}'s result from a truncated version of the stick-breaking construction (\eqn{\omega_{1d} = v_{1d}}; \eqn{\omega_{ld} = v_{ld}\prod_{r<l}(1-v_{dr})}, \eqn{l=2,\ldots,L_{d}}; \eqn{v_{d1},\ldots,v_{L_{d}-1}\sim} Beta \eqn{(1,\alpha_{d})}; \eqn{v_{Ld} = 1}, \eqn{\alpha_d \sim \Gamma(a_{\alpha_d},b_{\alpha_d})}), \eqn{\mathbf{\beta}_{ld}\sim N_{Q_d}(\mathbf{m}_{d},\mathbf{S}_{d})}, and \eqn{\sigma_{ld}^{-2}\sim\Gamma(a_{d},b_{d})}. It is further assumed that \eqn{\mathbf{m}_{d} \sim N_{Q_k}(\mathbf{m}_{0d},\mathbf{S}_{0d})} and \eqn{\mathbf{S}_{d}^{-1}\sim W(\nu,(\nu_k\Psi_d)^{-1})}. Here \eqn{W(\nu,(\nu\Psi)^{-1})} denotes a Wishart distribution with \eqn{\nu} degrees of freedom and expectation \eqn{\Psi^{-1}}, Here \eqn{\Gamma(a,b)} denotes a Gamma distribution with shape parameter \eqn{a} and rate parameter \eqn{b}, and \eqn{Q_d} denotes the dimension of the vector \eqn{\mathbf{z}_{di}}. It is worth mentioning that when \eqn{L_d=1}, the model for the conditional distribution of the test outcomes reduces to a normal regression model (where continuous covariates effects are modelled either parametrically or nonparametrically). For a detailed description, we refer to Inacio de Carvalho et al. (2013). 

The covariate-specific area under the curve is
\deqn{AUC(\mathbf{x})=\int_{0}^{1}ROC(p|\mathbf{x})dp.}
When the upper bound on the number of mixture components is 1, i.e., \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the covariate-specific AUC (binormal model), which is used in the package. In contrast, when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1}, the integral is computed numerically using Simpson's rule. With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is
\deqn{pAUC_{FPF}(u_1|\mathbf{x})=\int_0^{u_1} ROC(p|\mathbf{x})dp.} 
As for the AUC, when \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the \eqn{pAUC_{FPF}} (Hillis and Metz, 2012), and when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1} the integral is approximated numerically using Simpson's rule. The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1|\mathbf{x})/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as
\deqn{pAUC_{TPF}(u_2|\mathbf{x})=\int_{u_2}^{1}ROC_{TNF}(p|\mathbf{x})dp,}
where \eqn{ROC_{TNF}(p|\mathbf{x})} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as \eqn{ROC_{TNF}(p|\mathbf{x}) = F_{\bar{D}}\{F_{D}^{-1}(1-p|\mathbf{x})|\mathbf{x}\}.} Again, when \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the \eqn{pAUC_{TNF}} (Hillis and Metz, 2012), and when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1} the integral is approximated numerically using Simpson's rule. The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2|\mathbf{x})/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

It is worth referring that with respect to the computation of the DIC, when \eqn{L=1}, it is computed as in Spiegelhalter et al. (2002), and when \eqn{L>1}, DIC3 as described in Celeux et al. (2006) is computed. Also, for the computation of the conditional predictive ordinates (CPO) we follow the stable version proposed by Gelman et al. (2014).
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{newdata}{A data frame containing the values of the covariates at which the covariate-specific ROC curve (as well as the AUC, pAUC, dens and reg.fun, if required) was computed.}
\item{data}{The original supplied data argument.}
\item{missing.ind}{A logical value indicating whether for each pair of observations (test outcomes and covariates) missing values occur.}
\item{marker}{The name of the diagnostic test variable in the dataframe.}
\item{group}{The value of the argument \code{group} used in the call.}
\item{tag.h}{The value of the argument \code{tag.h} used in the call.}
\item{p}{Set of false positive fractions (FPF) at which the covariate-specific ROC curve has been estimated.}
\item{prior}{A list returning the hyperparameter values in the healthy and diseased populations.}
\item{ROC}{A list returning the \code{np} (length of the vector \code{p}) by npred predicted covariate-specific ROC curves (cROC) (posterior mean) and 95\% pointwise posterior credible bands.}
\item{AUC}{Estimated area under the covariate-specific ROC curve (posterior mean) and 95\% posterior credible band.}
\item{pAUC}{If computed, estimated partial area under the covariate-adjusted ROC curve (posterior mean) and 95\%  credible
band. Note that the returned values are normalised, so that the maximum value is one (see more on Details).}
\item{dens}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: \code{grid} (grid of test outcomes where the conditional densities were evaluated) and \code{dens} (MCMC realisations of the corresponding conditional densities).}
\item{reg.fun}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a data frame containing the predicted regression function (posterior mean) and 95\% credible band.}
\item{lpml}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: the log pseudo marginal likelihood (LPML) and the conditional predictive ordinates (CPO).}
\item{WAIC}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: widely applicable information criterion (WAIC) and associated complexity penalty (pW).}
\item{DIC}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: deviance information criterion (DIC) and associated complexity penalty (pD).}
\item{fit}{Named list of length two, with components \code{h} (healthy) and \code{d} (diseased). Each component is a list with the following information: (1) \code{formula}: the value of the argument \code{formula.h} or \code{formula.d} used in the call. (2) \code{mm}: information needed to construct the model matrix associated with single-weights linear dependent Dirichlet process mixture of normals model. (3) \code{beta}: array of dimension \code{nsave}x\code{L}x\code{Q} with the sampled regression coefficients. (4) \code{sd}: matrix of dimension \code{nsave}x\code{L} with the sampled variances. (4) \code{probs}: matrix of dimension \code{nsave}x\code{L} with the sampled components' weights. Here, \code{nsave} is the number of Gibbs sampler iterations saved, \code{L} is the maximum number of mixture components, and \code{Q} is the dimension of vector \eqn{\mathbf{z}_{d}}, \eqn{d \in \{D, \bar{D}\}}. (see also Details). (see also Details).}
\item{data_model}{A list with the data used in the fit: observed diagnostic test outcome and design matrices, separately for the healthy and diseased groups.}
}
\references{
Celeux, G., Forbes, F., Robert C. P., and Titerrington, D. M. (2006). Deviance information criteria for missing data models. Bayesian Analysis, \bold{1}, 651--674.

De Iorio, M., Johnson, W. O., Muller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. Biometrics, \bold{65}, 762--775.

Geisser, S. and Eddy, W.F. (1979) A Predictive Approach to Model Selection, Journal of the American Statistical Association, \bold{74}, 153--160.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and Rubin, D.B. (2014). Bayesian Data Analysis, 3rd ed. CRC Press: Boca Raton, FL.

Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Statistics and Computing, \bold{24}, 997--1010.

Hillis, S. L. and Metz, C.E. (2012). An Analytic Expression for the Binormal Partial Area under the ROC Curve. Academic Radiology, \bold{19}, 1491--1498.

Inacio de Carvalho, V., Jara, A., Hanson, T. E., and de Carvalho, M. (2013). Bayesian nonparametric ROC regression modeling. Bayesian Analysis, \bold{8}, 623--646.

Speigelhalter, D. J., Best, N. G., Carlin, B. P., and van der Linde, A. (2002). Bayesian measures of model comparison and fit. Journal of the Royal Statistical Society, Ser. B, \bold{64}, 583--639.

Watanabe, S. (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular Learning Theory. Journal of Machine Learning Research, \bold{11}, 3571--3594.
}
%\author{
%%  ~~who you are~~
%}
\note{
The input arguments \code{formula.h} and \code{formula.d} are similar to that used for the \code{\link{glm}} function, except that flexible specifications can be added by means of the function \code{f()}. For instance, specification \eqn{y \sim x1 + f(x2, K = 3)} would assume a linear effect of \code{x1} (if \code{x1} continuous) and the effect of \code{x2} would be modeled using B-splines basis functions. The argument \code{K = 3} indicates that \code{3} internal knots will be used, with the quantiles of \code{x2} used for their location. Categorical variables (factors) can be also incorporated, as well as interaction terms. For example, to include the factor-by-curve interaction between \code{age} and \code{gender} we need to specify, e.g., \eqn{y \sim gender + f(age, by = gender, K = c(3, 5))}. Note that, in this case, the number of knots can be different for each level of the factor. The order of the vector \code{K} of knots should match the levels of the factor.}

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
cROC_bnp <- cROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
               formula.d = l_marker1 ~ f(age, K = 0),
               group = "status", 
               tag.h = 0,
               data = newpsa,
               standardise = TRUE, 
               p = seq(0, 1, len = 101),
               compute.lpml = TRUE, 
               compute.WAIC = TRUE,
               compute.DIC = TRUE, 
               pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
               density = densitycontrol(compute = TRUE, grid.h = NA, grid.d = NA),
               mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1))
summary(cROC_bnp)

plot(cROC_bnp)
}
\dontshow{
agep <- seq(min(newpsa$age), max(newpsa$age), length = 5)
df.pred <- data.frame(age = agep)
cROC_bnp <- cROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
               formula.d = l_marker1 ~ f(age, K = 0),
               group = "status", 
               tag.h = 0,
               data = newpsa,
               newdata = df.pred,
               standardise = TRUE, 
               p = seq(0, 1, len = 101),
               compute.lpml = FALSE, 
               compute.WAIC = FALSE,
               compute.DIC = FALSE, 
               pauc = pauccontrol(compute = FALSE, value = 0.5, focus = "FPF"),
               density = densitycontrol(compute = FALSE, grid.h = NA, grid.d = NA),
               mcmc = mcmccontrol(nsave = 50, nburn = 10, nskip = 1))

summary(cROC_bnp)

plot(cROC_bnp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
