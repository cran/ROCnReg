\name{densitycontrol}
\alias{densitycontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
(Conditional) density estimates of test outcomes
}
\description{
This function is used to set various parameters controlling the estimation of the (conditional) density (densities) of test outcomes in both the healthy and diseased groups.
}
\usage{
densitycontrol(compute = FALSE, grid.h = NA, grid.d = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compute}{Logical value. If TRUE the (conditional) density (densities) of test outcomes in each group, healthy and diseased, are estimated.}
  \item{grid.h}{Grid of test outcomes in the healthy group where the (conditional) density (densities) estimates are to be evaluated. Value \code{NA} signals autoinitialization, with default a vector of length 200 in the range of test outcomes in the healthy group.}
  \item{grid.d}{Grid of test outcomes in the diseased group where the (conditional) density (densities) estimates are to be evaluated. Value \code{NA} signals autoinitialization, with default a vector of length 200 in the range of test outcomes in the diseased group.}
}
\details{
The value returned by this function is used as a control argument of the \code{\link{cROC.bnp}} and \code{\link{pooledROC.dpm}} functions.
}
\value{
A list with components for each of the possible arguments.
}

\seealso{
\code{\link{cROC.bnp}} and \code{\link{pooledROC.dpm}}
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

summary(cROC_bnp)

plot(cROC_bnp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
