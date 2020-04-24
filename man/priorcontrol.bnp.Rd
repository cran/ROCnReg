\name{priorcontrol.bnp}
\alias{priorcontrol.bnp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior information for the \code{AROC.bnp} and \code{cROC.bnp}
}
\description{
This function is used to set various parameters controlling the prior information to be used in the \code{AROC.bnp} and  \code{cROC.bnp} functions.
}
\usage{
priorcontrol.bnp(m0 = NA, S0 = NA, nu = NA, Psi = NA, a = 2, b = NA, 
	aalpha = 2, balpha = 2, L = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m0}{A numeric vector. Hyperparameter; mean vector of the (multivariate) normal prior distribution for the mean of the normal component of the centring distribution. \code{NA} signals autoinitialization, with defaults: a vector, of length \eqn{Q}, of zeros, if the data are standardised and the least squares estimates of the regression coefficients if the data are not standardised.}
  \item{S0}{A numeric matrix. Hyperparameter; covariance matrix of the (multivariate) normal prior distribution for the mean of the normal component of the centring distribution. \code{NA} signals autoinitialization, with defaults: 10\eqn{I_{Q\times Q}} if the data are standardised and \eqn{\mathbf{\hat{\Sigma}}} if the data are not standardised, where \eqn{\mathbf{\hat{\Sigma}}} is the asymptotic covariance matrix obtained of the regression coefficients obtained by fitting a linear model to the data.}
 \item{nu}{A numeric value. Hyperparameter; degrees of freedom of the Wishart prior distribution for the precision matrix of the the normal component of the centring distribution.\code{NA} signals autoinitialization, with default: \eqn{Q+2} where \eqn{Q} is the number of columns of the design matrix.}
  \item{Psi}{A numeric matrix. Hyperparameter; scale matrix of the Wishart distribution for the precision matrix of the the normal component of the centring distribution. \code{NA} signals autoinitialization, with defaults: \eqn{I_{Q\times Q}} if the data are standardised and to 30\eqn{\mathbf{\hat{\Sigma}}} if the data are not standardised, where \eqn{\mathbf{\hat{\Sigma}}} is the asymptotic covariance matrix obtained of the regression coefficients obtained by fitting a linear model to the data.}
 \item{a}{A numeric value. Hyperparameter; shape parameter of the gamma prior distribution for the precisions (inverse variances) of each component. The default is 2.}
  \item{b}{A numeric value. Hyperparameter; shape parameter of the gamma prior distribution for the precisions (inverse variances) of each component. \code{NA} signals autoinitialization, with defaults: 2 if the data are standardised and \eqn{\hat{\sigma}^2} if the data are not standardised}
  \item{aalpha}{A numeric value. Hyperparameter; shape parameter for the Gamma prior distribution of the precision parameter of the Dirichlet process prior. The default is 2.}
  \item{balpha}{A numeric value. Hyperparameter; rate parameter for the Gamma prior distribution of the precision parameter of the Dirichlet process prior. The default is 2.}
  \item{L}{A numeric value. Upper bound on the number of mixture components. Setting L = 1 corresponds to a normal model. The default is 10.}
}
%\details{
%TODO
%}
\value{
A list with components for each of the possible arguments.

}
%\references{
%% ~put references to the literature/web site here 
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AROC.bnp}} and  \code{\link{cROC.bnp}}
}
%\examples{
%# TODO
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
