waicnp <-
function(y, X, res, L, termsum = NULL) {
  n <- length(y)  
  p <- res$P
  beta <- res$Beta
  sigma2 <- res$Sigma2
  niter <- nrow(p)
    
  if(is.null(termsum)) {
    termsum <- inf_criteria(y, X, res)
  }
  logtermsum <- log(termsum)

  lpd <- sum(log(apply(exp(logtermsum),2,mean)))
  p2 <- sum(apply(logtermsum,2,var))
  waic <- -2*(lpd-p2)
  
  res <- list()
  res$pW <- p2
  res$WAIC <- waic
  res
}
