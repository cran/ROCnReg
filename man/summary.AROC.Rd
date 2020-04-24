\name{summary.AROC}
\alias{summary.AROC}
\title{
Summary method for \code{AROC} objects
}
\description{
Default summary method for objects fitted with \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()} functions.
}
\usage{
\method{summary}{AROC}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{AROC} as produced by \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()}.}
  \item{\dots}{Further arguments passed to or from other methods. Not yet implemented.}
}
\details{
%The information printed depends on the method. In all cases, the call to the function, the method, the area under the covariate-adjusted ROC curve (AAUC), the partial area under the covariate-adjusted ROC curve (if required) (AAUC), and the sample sizes are printed. For the semiparametric approach (\code{AROC.sp()}), the estimated coefficients (and 95\% confidence intervals, if required) of the model for the healthy population are printed. In addition, the function provides the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC). For the semiparametric Bayesian approach (\code{AROC.bsp()}), the posterior mean of the coefficients (and 95\% credible intervals) of the model for the healthy population are printed. For both Bayesian methods (\code{AROC.bsp()} and \code{AROC.bnp()}), and if required, the function provides the log pseudo marginal likelihood (LPML), the widely applicable information criterion (WAIC) and/or the deviance information criterion (DIC). For the kernel-based approach (\code{AROC.kernel()}), information regarding the selected bandwidth and the type of kernel estimator (for both regression and variance functions) is printed.
The information printed depends on the method. In all cases, the call to the function, the method, the area under the covariate-adjusted ROC curve (AAUC), the partial area under the covariate-adjusted ROC curve (if required) (AAUC), and the sample sizes are printed. For the semiparametric approach (\code{AROC.sp()}), the estimated coefficients (and 95\% confidence intervals, if required) of the model for the healthy population are printed. In addition, the function provides the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC). For the nonparametric Bayesian approach (\code{AROC.bnp()}), and if required, the function provides the log pseudo marginal likelihood (LPML), the widely applicable information criterion (WAIC) and/or the deviance information criterion (DIC). For the kernel-based approach (\code{AROC.kernel()}), information regarding the selected bandwidth and the type of kernel estimator (for both regression and variance functions) is printed.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AROC.bnp}}, \code{\link{AROC.sp}} or \code{\link{AROC.kernel}}.
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
m0 <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = TRUE, compute.WAIC = TRUE)

summary(m0)
}
\dontshow{
m0 <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
mcmc = list(nsave = 400, nburn = 100, nskip = 1))

summary(m0)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
