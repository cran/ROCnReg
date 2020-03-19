lpmlp <-
function(y, X, res, term = NULL) {
  n <- length(y)
  beta <- res$Beta
  sigma2 <- res$Sigma2
  niter <- length(sigma2)

  if(is.null(term)) {
	  term <- inf_criteria(y, X, res)
  }
  
  aux1 <- 1/term
  omegabari <- apply(aux1,2,mean)
  omegabari_1 <- sqrt(niter) * omegabari
  omegatilde <- matrix(0, nrow = niter, ncol = n)

  for(i in 1:n) {
    omegatilde[,i] <- pmin(aux1[,i], omegabari_1[i])
  }

  sum_omegatilde <- apply(omegatilde,2,sum)
  sum_term_omegatilde <- apply(term*omegatilde, 2, sum)
  cpo <- sum_term_omegatilde/sum_omegatilde

  lpml <- sum(log(cpo))

  res <- list()
  res$cpo <- cpo
  res$lpml <- lpml
  res 
}
