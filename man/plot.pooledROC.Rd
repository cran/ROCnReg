\name{plot.pooledROC}
\alias{plot.pooledROC}
\title{
Default pooledROC plotting
}
\description{
Takes a fitted \code{pooledROC} object produced by \code{pooledROC.BB}, \code{pooledROC.emp}, \code{pooledROC.kernel}, or \code{pooledROC.dpm} and plots the pooled ROC curve and associated area under the ROC curve (AUC).
}
\usage{
\method{plot}{pooledROC}(x, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{pooledROC} as produced by \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}} or \code{\link{pooledROC.dpm}}.}
  \item{main}{Character string with the overall title for the plot. If NULL, the default, the method used to estimate the pooled ROC curve is depicted.}
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
\code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}} or \code{\link{pooledROC.dpm}}}

\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)
\donttest{
m0_emp <- pooledROC.emp(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 500)

summary(m0_emp)

plot(m0_emp)
}
\dontshow{
m0_emp <- pooledROC.emp(marker = "l_marker1", group = "status",
tag.h = 0, data = newpsa, p = seq(0,1,l=101), B = 0)

summary(m0_emp)

plot(m0_emp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
