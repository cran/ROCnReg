\name{print.AROC}
\alias{print.AROC}
\title{
Print method for \code{AROC} objects
}
\description{
Default print method for objects fitted with \code{AROC.bnp()}, \code{AROC.sp()}, and \code{AROC.kernel()} functions.
}
\usage{
\method{print}{AROC}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{AROC} as produced by \code{AROC.bnp()}, \code{AROC.sp()}, or \code{AROC.kernel()}.}
  \item{\dots}{Further arguments passed to or from other methods. Not yet implemented.}
}
\details{
A short summary is printed including the area under the covariate-adjusted ROC curve (AAUC), and if required, the partial area under the covariate-adjusted ROC curve (pAAUC).
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
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = TRUE, compute.WAIC = TRUE)

AROC_bnp
}
\dontshow{
AROC_bnp <- AROC.bnp(formula.h = l_marker1 ~ f(age, K = 0),
group = "status", tag.h = 0, data = newpsa, standardise = TRUE,
p = seq(0,1,l=101), compute.lpml = FALSE, compute.WAIC = FALSE,
mcmc = list(nsave = 400, nburn = 100, nskip = 1))

AROC_bnp
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
