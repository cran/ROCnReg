\name{densitycontrol.aroc}
\alias{densitycontrol.aroc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional density estimates of test outcomes in the healthy population
}
\description{
This function is used to set various parameters controlling the estimation of the conditional densities of test outcomes in the healthy group.
}
\usage{
densitycontrol.aroc(compute = FALSE, grid.h = NA, newdata = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compute}{Logical value. If TRUE the conditional densities of test outcomes in the healthy group are estimated.}
  \item{grid.h}{Grid of test outcomes in the healthy group where the conditional density estimates are to be evaluated. Value \code{NA} signals autoinitialization, with default a vector of length 200 in the range of test outcomes in the healthy group.}
  \item{newdata}{Data frame containing the values of the covariates at which the conditional density estimates are computed.}
}
\details{
The value returned by this function is used as a control argument of the \code{\link{AROC.bnp}} function.
}
\value{
A list with components for each of the possible arguments.
}
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
\code{\link{AROC.bnp}} 
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)

# Covariate for prediction
agep <- seq(min(newpsa$age), max(newpsa$age), length = 5)
df.pred <- data.frame(age = agep)

\donttest{
AROC_bnp <- AROC.bnp(formula.h =  l_marker1 ~ f(age, K = 0),
                     group = "status", 
                     tag.h = 0,
                     data = newpsa,
                     standardise = TRUE,
                     p = seq(0, 1, len = 101),
                     compute.lpml = TRUE,
                     compute.WAIC = TRUE,
                     compute.DIC = TRUE,
                     pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
                     density = densitycontrol.aroc(compute = TRUE, grid.h = NA, newdata = df.pred),
                     mcmc = mcmccontrol(nsave = 500, nburn = 100, nskip = 1)
)
}
\dontshow{
AROC_bnp <- AROC.bnp(formula.h =  l_marker1 ~ f(age, K = 0),
                     group = "status", 
                     tag.h = 0,
                     data = newpsa,
                     standardise = TRUE,
                     p = seq(0, 1, len = 101),
                     compute.lpml = TRUE,
                     compute.WAIC = TRUE,
                     compute.DIC = TRUE,
                     pauc = pauccontrol(compute = TRUE, value = 0.5, focus = "FPF"),
                     density = densitycontrol.aroc(compute = TRUE, grid.h = NA, newdata = df.pred),
                     mcmc = mcmccontrol(nsave = 50, nburn = 10, nskip = 1)
)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
