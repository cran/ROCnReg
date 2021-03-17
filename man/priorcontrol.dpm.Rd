\name{priorcontrol.dpm}
\alias{priorcontrol.dpm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior information for the \code{pooledROC.dpm}
}
\description{
This function is used to set various parameters controlling the prior information to be used in the \code{pooledROC.dpm} function.
}
\usage{
priorcontrol.dpm(m0 = NA, S0 = NA, a = 2, b = NA, alpha = 1, L = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m0}{A numeric value. Hyperparameter; mean of the normal prior distribution for the means of each component. \code{NA} signals autoinitialization, with defaults: 0 if the data are standardised and \eqn{\bar{y}_d} (\eqn{d \in \{D, \bar{D}\}} if the data are not standardised.}
  \item{S0}{A numeric value. Hyperparameter; variance of the normal prior distribution for the means of each component. \code{NA} signals autoinitialization, with defaults: 10 if the data are standardised and 100*\eqn{\frac{s^2_d}{n_d}} (\eqn{d \in \{D, \bar{D}\}} if the data are not standardised, where \eqn{s_d} denotes the sample standard deviation.}
  \item{a}{A numeric value. Hyperparameter; shape parameter of the gamma prior distribution for the precisions (inverse variances) of each component. The default is 2.}
  \item{b}{A numeric value. Hyperparameter; rate parameter of the gamma prior distribution for the precisions (inverse variances) of each component. \code{NA} signals autoinitialization, with defaults: 0.5 if the data are standardised and \eqn{\frac{s^2_d}{2}} (\eqn{d \in \{D, \bar{D}\}} if the data are not standardised.}
  \item{alpha}{A numeric value. Precision parameter of the Dirichlet Process. The default is 1.}
  \item{L}{A numeric value. Upper bound on the number of mixture components. Setting L=1 corresponds to a normal model. The default is 10.}
}
%\details{
%TODO
%}
\value{
A list with components for each of the possible arguments.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pooledROC.dpm}}
}
%\examples{
%# TODO
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
