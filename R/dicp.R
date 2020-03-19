dicp <-
function(y, X, res, term = NULL) {
  n <- length(y)
  beta <- res$Beta
  sigma2 <- res$Sigma2
  niter <- length(sigma2)

  if(is.null(term)) {
	  term <- inf_criteria(y, X, res)
  }

  logterm <- log(term)
  logtermsum <- apply(logterm, 1, sum)

  deviance <- -2*mean(logtermsum)
  
  betahat <- apply(beta,2,mean)
  sigma2hat <- mean(sigma2)
  
  d <- -2*sum(dnorm(y, X%*%betahat, sqrt(sigma2hat), log = TRUE))
  pd <- deviance - d
  DIC <- deviance + pd 
  
  res <- list()
  res$pD <- pd
  res$DIC <- DIC
  res
}
