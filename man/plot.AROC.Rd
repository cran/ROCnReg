\name{plot.AROC}
\alias{plot.AROC}
\title{
Default AROC plotting
}
\description{
Takes a fitted \code{AROC} object produced by \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()} and plots the covariate-adjusted ROC curve (AROC) and associated area under the AROC (AAUC).
}
\usage{
\method{plot}{AROC}(x, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{AROC} as produced by \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()}.}
   \item{main}{Character string with the overall title for the plot. If NULL, the default, the method used to estimate the AROC curve is depicted.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
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
\code{\link{AROC.bnp}}, \code{\link{AROC.sp}} or \code{\link{AROC.kernel}}}
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
p = seq(0,1,l=101), compute.lpml = TRUE, compute.WAIC = TRUE,
compute.DIC = TRUE)

plot(AROC_bnp)
}
\dontshow{
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
compute.DIC = FALSE,
mcmc = list(nsave = 400, nburn = 100, nskip = 1))

plot(AROC_bnp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
