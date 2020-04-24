\name{mcmccontrol}
\alias{mcmccontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Markov chain Monte Carlo (MCMC) parameters
}
\description{
This function is used to set various parameters controlling the Markov chain Monte Carlo (MCMC) parameters.
}
\usage{
mcmccontrol(nsave = 8000, nburn = 2000, nskip = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsave}{An integer giving the total number of scans to be saved (does not include the burn-in and thinning iterations). }
  \item{nburn}{An integer giving the number of burn-in scans.}
  \item{nskip}{An integer giving the thinning interval,}
}
\details{
The value returned by this function is used as a control argument of the \code{\link{AROC.bnp}}, \code{\link{cROC.bnp}}, and \code{\link{pooledROC.dpm}} functions.
}
\value{
A list with components for each of the possible arguments.
}
%\references{
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pooledROC.dpm}}, \code{\link{AROC.bnp}} and \code{\link{cROC.bnp}}
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
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
