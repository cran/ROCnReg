dic <-
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
  deviance <- -2* mean(apply(logtermsum, 1, sum)) 
  pd <- deviance + 2 * sum(log(apply(exp(logtermsum),2,mean)))
  DIC <- deviance + pd
  
  res <- list()
  res$pD <- pd
  res$DIC <- DIC
  res
}
