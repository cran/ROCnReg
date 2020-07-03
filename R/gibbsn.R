gibbsn <-
function(y, prior, mcmc, standardise = FALSE) {
  n <- length(y)
  
  yt <- y
  if(standardise) {
    yt <- (y-mean(y))/sd(y)
  }  
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nsave*nskip+nburn
  
  m0 <- prior$m0
  S0 <- prior$S0
  a <- prior$a 
  b <- prior$b
  
  Mu <- Mu1 <- Sigma2 <- Sigma21 <- numeric(nsim)
  Mu[1] <- mean(yt)
  Sigma2[1] <- var(yt)
  
  for(i in 2:nsim) {
    meanmu <- ((n*mean(yt)*S0+m0*Sigma2[i-1])/(n*S0+Sigma2[i-1]))
    varmu <- (Sigma2[i-1]*S0)/(n*S0+Sigma2[i-1])
    Mu[i] <- Mu1[i] <- rnorm(1,meanmu,sqrt(varmu))
    if(standardise) {
      Mu1[i] <- sd(y)*Mu[i]+mean(y)
    }
    
    a1 <- a+n/2
    b1 <- b+0.5*sum((yt-Mu[i])^2)
    Sigma2[i] <- Sigma21[i] <- 1/rgamma(1,a1,b1)
    if(standardise) {
      Sigma21[i] <- var(y)*Sigma2[i]
    }
  }
  
  res <- list()
  res$Mu <- Mu1[seq(nburn+1, nsim, by = nskip)]
  res$Sigma2 <- Sigma21[seq(nburn+1, nsim, by = nskip)]
  return(res)
}
