lpml <-
function(y, X, res, L, termsum = NULL) {
  n <- length(y)
  p <- res$P
  beta <- res$Beta
  sigma2 <- res$Sigma2
  niter <- nrow(p)
 
  if(is.null(termsum)) {
    termsum <- inf_criteria(y, X, res)
  }
  
  aux <- 1/termsum
  omegabari <- apply(aux, 2, mean)
  omegabari_1 <- sqrt(niter) * omegabari
  omegatilde <- matrix(0, nrow = niter, ncol = n)
  
  for(i in 1:n) {
    omegatilde[,i] <- pmin(aux[,i], omegabari_1[i])  
  }

  sum_omegatilde <- apply(omegatilde,2,sum)
  sum_term_omegatilde <- apply(termsum*omegatilde, 2, sum)
  cpo <- sum_term_omegatilde/sum_omegatilde
 
  lpml <- sum(log(cpo))
  
  res <- list()
  res$cpo <- cpo
  res$lpml <- lpml
  res 
}
