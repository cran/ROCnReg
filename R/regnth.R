regnth <-
function(y, X, prior, mcmc, standardise = TRUE) {
    yt <- y
    if(standardise == TRUE){
        yt <- (y-mean(y))/sd(y)
    }
    n <- length(y)
    k <- ncol(X)
    
    m0 <- prior$m0
    S0 <- prior$S0
    nu0 <- prior$nu
    psi0 <- prior$Psi
    a <- prior$a
    b <- prior$b
    
    nburn <- mcmc$nburn
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nsim <- nburn + nskip * nsave
    
    beta <- beta1 <- matrix(0, nrow = nsim, ncol = k, dimnames = list(1:nsim, colnames(X)))
    sigma2 <- sigma21 <- numeric(nsim)
    mu0 <- matrix(0, nrow = nsim, ncol = k)
    V0inv <- array(0, c(nsim,k,k))
    
    bols <- ols.function(X, yt, vcov = FALSE)$coeff #solve(t(X)%*%X)%*%(t(X)%*%yt)
    e <- yt-X%*%bols
    s2 <- crossprod(e)/(n-k)
    beta[1,] <- bols
    sigma2[1] <- s2
    if(standardise == TRUE){
        beta1[1,1] <- sd(y)*beta[1,1] + mean(y)
        if(k > 1) {
            beta1[1,2:k] <- sd(y)*beta[1,2:k]
        }
        sigma21 <- var(y)*sigma2[1]
    }
    mu0[1,] <- rep(0,k)
    V0inv[1,,] <- solve(100*diag(k))
    
    for(i in 2:nsim) {
        V1 <- solve(V0inv[i-1,,]+(1/sigma2[i-1])*crossprod(X))
        mu1 <- V1%*%(V0inv[i-1,,]%*%mu0[i-1,]+(1/sigma2[i-1])*crossprod(X,yt))
        beta1[i,] <- beta[i,] <- mvrnorm(1, mu = mu1, Sigma = V1) #rmvn(1, mu = mu1, sigma = V1, isChol = TRUE)
        if(standardise == TRUE){
            beta1[i,1] <- sd(y)*beta[i,1] + mean(y)
            if(k > 1) {
                beta1[i,2:k] <- sd(y)*beta[i,2:k]
            }
        }
        
        a1 <- a + (n/2)
        b1 <- b + 0.5*crossprod(yt-X%*%beta[i,])
        sigma21[i] <- sigma2[i] <- 1/rgamma(1, shape = a1, rate = b1)
        if(standardise == TRUE){
            sigma21[i] <- var(y)*sigma2[i]
        }
        Vaux <- solve(solve(S0) + V0inv[i-1,,])
        mu0[i,] <- mvrnorm(1, mu = Vaux%*%(V0inv[i-1,,]%*%beta[i,]+solve(S0)%*%m0), Sigma = Vaux) 
                   #mvrnorm(1, mu = Vaux%*%(V0inv[i-1,,]%*%t(t(beta[i,]))+solve(S0)%*%m0), Sigma = Vaux) 
                   #rmvn(1, mu = Vaux%*%(V0inv[i-1,,]%*%t(t(beta[i,]))+solve(S0)%*%m0), sigma = Vaux)
        
        V0inv[i,,] <- rWishart(1, df = nu0+1, solve(nu0*psi0 + as.numeric(crossprod(beta[i,]-mu0[i,]))))
        #V0inv[i,,] <- rWishart(1, df = nu0+1,solve(nu0*psi0+as.numeric((beta[i,]-mu0[i,])%*%t(t(beta[i,]-mu0[i,])))))
        
    }
    res <- list()
    res$Beta <- beta1[seq(nburn+1, nsim, by = nskip),,drop = FALSE]
    res$Sigma2 <- sigma21[seq(nburn+1, nsim, by = nskip)]
    res
}
