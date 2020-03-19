waicp <-
function(y, X, res, term = NULL) {
  n <- length(y)
  beta <- res$Beta
  sigma2 <- res$Sigma2
  niter <- length(sigma2)

  if(is.null(term)) {
	  term <- inf_criteria(y, X, res)
  }

  logterm <- log(term)
  lpd <- sum(log(apply(exp(logterm),2,mean)))
  p2 <- sum(apply(logterm,2,var))
  waic <- -2*(lpd-p2)

  res <- list()
  res$pW <- p2
  res$WAIC <- waic
  res
}
