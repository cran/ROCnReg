\name{pauccontrol}
\alias{pauccontrol}
\title{
Partial area under the covariate-adjusted/covariate-specific/pooled ROC curve
}
\description{
Used to set various parameters controlling the estimation of the partial area under the covariate-adjusted/covariate-specific/pooled ROC curve .
}
\usage{
pauccontrol(compute = FALSE, focus = c("FPF", "TPF"), value = 1)
}
\arguments{
  \item{compute}{Logical value. If TRUE the partial area under the covariate-adjusted/covariate-specific/pooled ROC curve is estimated.}
  \item{focus}{Whether computation should be done over a restricted range of false positive fractions (FPF) or a restricted range of true positive fractions (TPF).}
  \item{value}{Numeric value. Pre-specified upper bound for the FPF (if focus = "FPF") or lower bound for the TPF (if focus = "TPF").}
}
\details{The value returned by this function is used as a control argument of the \code{\link{AROC.bnp}}, \code{\link{AROC.sp}}, \code{\link{AROC.kernel}}, \code{\link{cROC.bnp}}, \code{\link{cROC.sp}}, \code{\link{cROC.kernel}} functions.
}
\value{A list with components for each of the possible arguments.}

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
cROC_sp_normal <- cROC.sp(formula.h = l_marker1 ~ age,
                          formula.d = l_marker1 ~ age,
                          group = "status", 
                          tag.h = 0,
                          data = newpsa,
                          newdata = df.pred,
                          est.cdf = "normal",
                          pauc = list(compute = TRUE, value = 0.5, focus = "FPF"),
                          p = seq(0, 1, l = 101), 
                          B = 10)
}
\dontshow{
cROC_sp_normal <- cROC.sp(formula.h = l_marker1 ~ age,
                          formula.d = l_marker1 ~ age,
                          group = "status", 
                          tag.h = 0,
                          data = newpsa,
                          newdata = df.pred,
                          est.cdf = "normal",
                          pauc = list(compute = TRUE, value = 0.5, focus = "FPF"),
                          p = seq(0, 1, l = 101), 
                          B = 10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
