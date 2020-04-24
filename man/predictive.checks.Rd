\name{predictive.checks}
\alias{predictive.checks}
\title{
Posterior predictive checks.
}
\description{
Implements posterior predictive checks for objects of \code{AROC} or \code{cROC} as produced by \code{AROC.bnp} or \code{cROC.bnp}.
}
\usage{
predictive.checks(object, 
	statistics = c("min", "max", "kurtosis", "skewness"), 
	ndensity = 512, devnew = TRUE)
}
\arguments{
 	\item{object}{An object of class \code{AROC} or \code{cROC} as produced by \code{AROC.bnp} or \code{cROC.bnp}.}
  	\item{statistics}{Character vector. Statistics to be used for the posterior predictive checking. By default, "min", "max", "kurtosis" and "skewness"}
  	\item{ndensity}{An integer giving the number of equally spaced points at which the density of the test outcomes and of the simulated datasets from the posterior predictive distribution (see more on Details) is to be estimated (for more details see the help of the function \code{density} in the \code{stats} package).}
  	\item{devnew}{A logical value. If TRUE, each plot is depicted in a new graphic device.}
}
\details{
Compares a selected test statistic computed based on the observed diagnostic test outcome (either in the nondiseased group, \code{AROC} object, or in both the nondiseased and diseased groups, \code{cROC} object) against the same test statistics computed based on simulated data from the posterior predictive distribution of the diagnostic test outcome obtained using a single-weights linear dependent Dirichlet process mixture of normals model.
The following graphics are depicted: (1) histograms of the desired statistics computed from simulated datasets (nsave of them) from the posterior predictive distribution of the diagnostic test outcome. In these plots, the estimated statistics from the observed diagnostic test outcome are also depicted. (2) Kernel density estimates computed computed from simulated datasets (nsave of them) from the posterior predictive distribution of the diagnostic test outcome. In these plots, the kernel density estimate of the observed diagnostic test outcome is also depicted. In the case of an object of class \code{AROC}, the abovementioned graphics are depicted only for the diagnostic test outcome in the nondiseased group. However, for an object of class \code{cROC}, the graphics are depicted, separately, for both the nondiseased and diseased groups. For a detailed discussion about predictive checks, see Gabry et al. (2019).
}
\value{
As a result, the function provides a list with the following components:
\item{yrep}{List of matrices associated with the diseased (d) and nondiseased (h) groups. Each column of the matrix (there are nsave of them) corresponds to a dataset generated from the posterior predictive distribution of the diagnostic test outcomes. For \code{AROC} objects, the list only contains results for the nondiseased group.}
\item{y}{List of numeric vectors associated with the diseased (d) and nondiseased (h) groups. The vector contains the observed diagnostic test outcomes. For \code{AROC} objects, the list only contains results for the nondiseased group.}
}
\references{
Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society, Series A, 182, 1--14.

Inacio de Carvalho, V., and Rodriguez-Alvarez, M. X. (2018). Bayesian nonparametric inference for the covariate-adjusted ROC curve. arXiv preprint arXiv:1806.00473.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AROC.bnp}}, \code{\link{cROC.bnp}}
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
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = TRUE, compute.WAIC = TRUE)

predictive.checks(AROC_bnp, statistics = "skewness")
}
\dontshow{
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
mcmc = list(nsave = 400, nburn = 100, nskip = 1))

predictive.checks(AROC_bnp, statistics = "skewness")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
