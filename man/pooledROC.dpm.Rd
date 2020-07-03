\name{pooledROC.dpm}
\alias{pooledROC.dpm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Bayesian inference of the pooled ROC curve 
}
\description{
This function estimates the pooled ROC curve using a Dirichlet process mixture of normals model as proposed by Erkanli et al. (2006).
}
\usage{
pooledROC.dpm(marker, group, tag.h, data, 
  standardise = TRUE, p = seq(0, 1, l = 101),
  compute.lpml = FALSE, compute.WAIC = FALSE, compute.DIC = FALSE, 
  pauc = pauccontrol(), density = densitycontrol(), 
  prior.h = priorcontrol.dpm(), prior.d = priorcontrol.dpm(),
  mcmc = mcmccontrol(),
  parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{A character string with the name of the diagnostic test variable.}
  \item{group}{A character string with the name of the variable that distinguishes healthy from diseased individuals.}
  \item{tag.h}{The value codifying healthy individuals in the variable \code{group}.}
  \item{data}{Data frame representing the data and containing all needed variables.}
  \item{standardise}{A logical value. If TRUE the test outcomes are standardised (so that the resulting test outcomes have mean zero and standard deviation of one). The default is TRUE.}
  \item{p}{Set of false positive fractions (FPF) at which to estimate the pooled ROC curve. This set is also used to compute the area under the ROC curve (AUC) using Simpson's rule. Thus, the length of the set should be an odd number, and it should be rich enough for an accurate estimation.}
  \item{compute.lpml}{A logical value. If TRUE, the log pseudo marginal likelihood (LPML, Geisser and Eddy, 1979) and the conditional predictive ordinates (CPO) are computed.}
  \item{compute.WAIC}{A logical value. If TRUE, the widely applicable information criterion (WAIC, Gelman et al., 2014; Watanabe, 2010) is computed.}
  \item{compute.DIC}{A logical value. If TRUE, the deviance information criterion is computed.}
   \item{pauc}{A list of control values to replace the default values returned by the function \code{\link{pauccontrol}}. This argument is used to indicate whether the partial area under the pooled ROC curve should be computed, and in case it is computed, whether the focus should be placed on restricted false positive fractions (FPFs) or on restricted true positive fractions (TPFs), and the upper bound for the FPF (if focus is FPF) or the lower bound for the TPF (if focus is TPF).}
   \item{density}{A list of control values to replace the default values returned by the function \code{\link{densitycontrol}}. This argument is used to indicate whether the densities of the marker in the healthy and diseased population should be computed, and in case it is to be computed, at which grid of test outcomes in each of the populations.}
  \item{prior.h}{Hyparameter specification for the healthy population. A list of control values to replace the default values returned by the function \code{\link{priorcontrol.dpm}}. See \code{\link{priorcontrol.dpm}} for details.}
  \item{prior.d}{Hyparameter specification for the diseased population. A list of control values to replace the default values returned by the function \code{\link{priorcontrol.dpm}}. See \code{\link{priorcontrol.dpm}} for details.}
  \item{mcmc}{A list of control values to replace the default values returned by the function \code{\link{mcmccontrol}}. See \code{\link{mcmccontrol}} for details.}
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
The method implemented in this function estimates \eqn{F_{D}(\cdot)} and \eqn{F_{\bar{D}}(\cdot)} by means of a Dirichlet process mixture of normals model. More precisely, and letting \eqn{\{y_{\bar{D}i}\}_{i=1}^{n_{\bar{D}}}}  and \eqn{\{y_{Dj}\}_{j=1}^{n_{D}}} be two independent random samples from the nondiseased and diseased populations, respectively, the model postulated for the distribution function is as follows
\deqn{F_{\bar{D}}(y_{\bar{D}i}) = \sum_{l=1}^{L_{\bar{D}}}\omega_{l\bar{D}}\Phi(y_{\bar{D}i}\mid\mu_{l\bar{D}},\sigma_{l\bar{D}}^2),}
\deqn{F_{D}(y_{Dj}) = \sum_{l=1}^{L_{D}}\omega_{lD}\Phi(y_{Dj}\mid\mu_{lD},\sigma_{lD}^2),}
where \eqn{L_{d}} is pre-specified is a pre-specified upper bound on the number of mixture components (\eqn{d \in \{D, \bar{D}\}}). The \eqn{\omega_{ld}}'s result from a truncated version of the stick-breaking construction (\eqn{\omega_{1d} = v_{1d}}; \eqn{\omega_{ld} = v_{ld}\prod_{r<l}(1-v_{dr})}, \eqn{l=2,\ldots,L_{d}}; \eqn{v_{d1},\ldots,v_{L_{d}-1}\sim} Beta \eqn{(1,\alpha_{d})}; \eqn{v_{Ld} = 1}, \eqn{\alpha_d \sim \Gamma(a_{\alpha_d},b_{\alpha_d})}), \eqn{\beta_{ld}\sim N(m_{0d},S_{0d})}, and \eqn{\sigma_{ld}^{-2}\sim\Gamma(a_{d},b_{d})}.

The area under the curve is 
\deqn{AUC=\int_{0}^{1}ROC(p)dp.} 
When the upper bound on the number of mixture components is 1, i.e., \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the AUC (binormal model), which is used in the package. In contrast, when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1}, the AUC is computed using results presented in Erkanli et al. (2006). With regard to the partial area under the curve, when \code{focus = "FPF"} and assuming an upper bound \eqn{u_1} for the FPF, what it is computed is \deqn{pAUC_{FPF}(u_1)=\int_0^{u_1} ROC(p)dp.} 
As for the AUC, when \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the \eqn{pAUC_{FPF}} (Hillis and Metz, 2012), and when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1} the integral is approximated numerically using Simpson's rule. The returned value is the normalised pAUC, \eqn{pAUC_{FPF}(u_1)/u_1} so that it ranges from \eqn{u_1/2} (useless test) to 1 (perfect marker). Conversely, when \code{focus = "TPF"}, and assuming a lower bound for the TPF of \eqn{u_2}, the partial area corresponding to TPFs lying in the interval \eqn{(u_2,1)} is computed as \deqn{pAUC_{TPF}(u_2)=\int_{u_2}^{1}ROC_{TNF}(p)dp,} where \eqn{ROC_{TNF}(p)} is a \eqn{270^\circ} rotation of the ROC curve, and it can be expressed as \eqn{ROC_{TNF}(p) = F_{\bar{D}}\{F_{D}^{-1}(1-p)\}.} Again, when \eqn{L_d = 1} (\eqn{d \in \{D, \bar{D}\}}), there is a closed-form expression for the \eqn{pAUC_{TNF}} (Hillis and Metz, 2012), and when \eqn{L_D > 1} or \eqn{L_{\bar{D}} > 1} the integral is approximated numerically using Simpson's rule. The returned value is the normalised pAUC, \eqn{pAUC_{TPF}(u_2)/(1-u_2)}, so that it ranges from \eqn{(1-u_2)/2} (useless test) to 1 (perfect test).

It is worth referring that with respect to the computation of the DIC, when \eqn{L=1}, it is computed as in Spiegelhalter et al. (2002), and when \eqn{L>1}, DIC3 as described in Celeux et al. (2006) is computed. Also, for the computation of the conditional predictive ordinates (CPO) we follow the stable version proposed by Gelman et al. (2014).
}
\value{
As a result, the function provides a list with the following components:
\item{call}{The matched call.}
\item{marker}{A list with the diagnostic test outcomes in the healthy (h) and diseased (d) groups.}
\item{missing.ind}{A logical value indicating whether missing values occur.}
\item{p}{Set of false positive fractions (FPF) at which the pooled ROC curve has been estimated.}
\item{prior}{A list returning the hyperparameter values in the healthy and diseased populations.}
\item{ROC}{Estimated pooled ROC curve, and corresponding 95\% pointwise credible band.}
\item{AUC}{Estimated pooled AUC, and corresponding 95\% credible interval.}
\item{pAUC}{If computed, estimated partial area under the pooled ROC curve (posterior mean) and 95\% pointwise credible band. Note that the returned values are normalised, so that the maximum value is one (see more on Details).}
\item{dens}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: \code{grid} (grid of test outcomes where the densities are evaluated) and \code{dens} (MCMC realisations of the corresponding densities).}
\item{lpml}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: the log pseudo marginal likelihood (LPML) and the conditional predictive ordinates (CPO).}
\item{WAIC}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: widely applicable information criterion (WAIC) and associated complexity penalty (pW).}
\item{DIC}{If computed, named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with two components: deviance information criterion (DIC) and associated complexity penalty (pD).}
\item{fit}{Named list of length two, with components 'h' (healthy) and 'd' (diseased). Each component is a list with the following information: (1)\code{P}: matrix of dimension nsavexL with the sampled components' weights; (2) \code{Mu}: matrix of dimension nsavexL with the sampled means; and (3) \code{Sigma2}: matrix of dimension nsavexL with the sampled variances. Here, nsave is the number of Gibbs sampler iterates saved, and L is the upper bound on the number of mixture components.}
}
\references{
Erkanli, A., Sung M., Jane Costello, E., and Angold, A. (2006). Bayesian semi-parametric ROC analysis. Statistics in Medicine, \bold{25}, 3905--3928.

Geisser, S. and Eddy, W.F. (1979) A Predictive Approach to Model Selection, Journal of the American Statistical Association, \bold{74}, 153--160.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and Rubin, D.B. (2014). Bayesian Data Analysis, 3rd ed. CRC Press: Boca Raton, FL.

Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. Statistics and Computing, \bold{24}, 997--1010.

Hillis, S. L. and Metz, C.E. (2012). An Analytic Expression for the Binormal Partial Area under the ROC Curve. Academic Radiology, \bold{19}, 1491--1498.

Watanabe, S. (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular Learning Theory. Journal of Machine Learning Research, \bold{11}, 3571--3594.
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
m0_dpm <- pooledROC.dpm(marker = "l_marker1", group = "status",
            tag.h = 0, data = newpsa, standardise = TRUE, 
            p = seq(0,1,l=101), compute.WAIC = TRUE, compute.lpml = TRUE, 
            compute.DIC = TRUE, 
            prior.h = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            prior.d = priorcontrol.dpm(m0 = 0, S0 = 100, a = 2, b = 2, aalpha = 2, 
            balpha = 2, L =10),
            mcmc = mcmccontrol(nsave = 400, nburn = 100, nskip = 1))

summary(m0_dpm)

plot(m0_dpm)
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
            mcmc = mcmccontrol(nsave = 400, nburn = 100, nskip = 1))

summary(m0_dpm)

plot(m0_dpm)
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
