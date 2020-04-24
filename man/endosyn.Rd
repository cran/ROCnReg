\name{endosyn}
\alias{endosyn}
\docType{data}
\title{
Simulated endocrine data.
}
\description{
The \code{endosyn} data set was simulated based on the data analysed in Rodriguez-Alvarez et al. (2011a,b) and Inacio de Carvalho and Rodriguez-Alvarez (2018); and presented in Botana et al. (2007) and Tome et al. (2008). The aim of these studies was to use the body mass index (BMI) to detect patients having a higher risk of cardiovascular problems, ascertaining the possible effect of age and gender on the accuracy of this measure.
}
\usage{data("endosyn")}
\format{
  A data frame with 2840 observations on the following 4 variables.
  \describe{
  	\item{\code{gender}}{patient's gender. Factor with \code{Men} and \code{Women} levels.}
  	\item{\code{age}}{patient's age.}
    \item{\code{cvd_idf}}{true disease status (presence/absence of two of more cardiovascular risk factors according to the International Diabetes Federation). Numerical vector (0 = absence, 1 = presence).}
    \item{\code{bmi}}{patient's body mass index.}
  }
}

\source{
Botana, M.A., Mato, J.A., Cadarso-Suarez, C., Tome, M.A., Perez-Fernandez, R., Fernandez-Mario, A., Rego-Iraeta, A., Solache, I. (2007). Overweight, obesity and central obesity prevalences in the region of Galicia in Northwest Spain. Obesity and Metabolism, 3, 106--115.

Tome, M.A., Botana, M.A., Cadarso-Suarez, C., Rego-Iraeta, A., Fernandez-Mario, A., Mato, J.A, Solache, I., Perez-Fernandez, R. (2008). Prevalence of metabolic syndrome in Galicia (NW Spain) on four alternative definitions and association with insulin resistance. Journal of Endocrinological Investigation, 32, 505--511.

}
\references{
Inacio de Carvalho, V., and Rodriguez-Alvarez, M. X. (2018). Bayesian nonparametric inference for the covariate-adjusted ROC curve. arXiv preprint arXiv:1806.00473.

Rodriguez-Alvarez, M.X., Roca-Pardinas, J. and Cadarso-Suarez, C. (2011a). ROC curve and covariates: extending induced methodology to the non-parametric framework. Statistics and Computing, 21(4), 483--49.

Rodriguez- Alvarez, M.X., Roca-Pardinas, J. and Cadarso-Suarez, C. (2011b). A new flexible direct ROC regression model - Application to the detection of cardiovascular risk factors by anthropometric measures. Computational Statistics and Data Analysis, 55(12), 3257--3270.
}
\examples{
data(endosyn)
summary(endosyn)
}
\keyword{datasets}