\name{summary.pooledROC}
\alias{summary.pooledROC}
\title{
Summary method for \code{pooledROC} objects
}
\description{
Default summary method for objects fitted with \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}}, or \code{\link{pooledROC.dpm}} functions.
}
\usage{
\method{summary}{pooledROC}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{pooledROC} as produced by \code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}}, or \code{\link{pooledROC.dpm}}.}
  \item{\dots}{Further arguments passed to or from other methods. Not yet implemented.}
}
\details{
A short summary is printed including the call to the function, the method, samples sizes, the area under the pooled ROC curve (AUC), and if required, the partial area under the pooled ROC curve. For the nonparametric Bayesian approach (\code{pooledROC.dpm()}), and if required, the function provides the log pseudo marginal likelihood (LPML), the widely applicable information criterion (WAIC) and/or the deviance information criterion (DIC). For the kernel-based approach (\code{pooledROC.dpm()}), information regarding the selected bandwidths and the density bandwidth selection method is presented. 
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
\code{\link{pooledROC.BB}}, \code{\link{pooledROC.emp}}, \code{\link{pooledROC.kernel}} or \code{\link{pooledROC.dpm}}.
}
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
