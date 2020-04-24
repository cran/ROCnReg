\name{plot.cROC}
\alias{plot.cROC}
\title{
Default cROC plotting
}
\description{
Takes a fitted \code{cROC} object produced by \code{cROC.bnp()}, \code{cROC.sp()}, or \code{cROC.kernel()} and plots the covariate-specific ROC curve (cROC) and associated area under the cROC curve. The suitable type of graphic is chosen according to the number and nature of the covariates.
}
\usage{
\method{plot}{cROC}(x, ask = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{cROC} as produced by \code{cROC.bnp()}, \code{cROC.sp()}, or \code{cROC.kernel()}.}
  \item{ask}{A logical value. If TRUE, the default, the user is asked for confirmation, before a new figure is drawn.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%\details{
%TODO
%}
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
\code{\link{cROC.bnp}}, \code{\link{cROC.sp}} or \code{\link{cROC.kernel}}
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

plot(cROC_bnp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
