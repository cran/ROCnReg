\name{summary.cROC}
\alias{summary.cROC}
\title{
Summary method for \code{cROC} objects
}
\description{
Default summary method for objects fitted with \code{cROC.bnp()}, \code{cROC.sp()}, or \code{cROC.kernel()} functions.
}
\usage{
\method{summary}{cROC}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{cROC} as produced by \code{cROC.bnp()}, \code{cROC.sp()}, or \code{cROC.kernel()}.}
  \item{\dots}{Further arguments passed to or from other methods. Not yet implemented.}
}
\details{
The information printed depends on the method. In all cases, the call to the function, the method, and the sample sizes are printed. For the semiparametric approach (\code{cROC.sp()}), the estimated coefficients (and 95\% confidence intervals, if required) of the model for the healthy population, the diseased population and the conditional ROC curve, are printed. In addition, the function provides the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC). For the nonparametric Bayesian approach (\code{cROC.bnp()}), and if required, the function provides the log pseudo marginal likelihood (LPML), the widely applicable information criterion (WAIC) and/or the deviance information criterion (DIC) (for both healthy and diseased populations). For the kernel-based approach (\code{cROC.kernel()}), information regarding the selected bandwidths and the type of kernel estimator(for both healthy and diseased populations and for both regression and variance functions) is printed.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Rodriguez-Alvarez, M.X., Tahoces, P.G., Cadarso-Suarez, C. and Lado, M.J. (2011). Comparative study of ROC regression techniques. Applications for the computer-aided diagnostic system in breast cancer detection. \emph{Computational Statistics and Data Analysis}, \bold{55}, 888--902.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cROC.bnp}}, \code{\link{cROC.sp}} or \code{\link{cROC.kernel}}.
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
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
