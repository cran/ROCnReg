\name{cROCData}
\alias{cROCData}
\title{
Selects an adequate set of points from a data set for obtaining predictions.
}
\description{
Selects an adequate set of points from a data set for obtaining predictions 
}
\usage{
cROCData(data, names.cov, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data set from which the new set of covariate values is obtained.}
  \item{names.cov}{Character vector with the names of the covariates to be included in the new data set.}
  \item{group}{A character string with the name of the variable in the original data set that distinguishes healthy from diseased individuals.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A data frame containing selected values of all needed covariates. For those that are continuous, 30 different values are selected.
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
\code{\link{AROC.bnp}}, \code{\link{cROC.bnp}}, \code{\link{cROC.sp}}, \code{\link{cROC.kernel}}, \code{\link{compute.threshold.cROC}} or \code{\link{compute.threshold.AROC}}.
}
\examples{
library(ROCnReg)
data(psa)
# Select the last measurement
newpsa <- psa[!duplicated(psa$id, fromLast = TRUE),]

# Log-transform the biomarker
newpsa$l_marker1 <- log(newpsa$marker1)

newdf <- cROCData(newpsa, "age", "status")

summary(newdf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
