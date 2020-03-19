dpm <-
function(y, prior, mcmc, standardise = FALSE) {
    n <- length(y)
    yt <- y
    if(standardise) {
        yt <- (y-mean(y))/sd(y)
    }
    
    m0 <- prior$m0
    S0 <- prior$S0
    a <- prior$a
    b <- prior$b
    aalpha <- prior$aalpha
    balpha <- prior$balpha
    L <- prior$L
    
    nburn <- mcmc$nburn
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nsim <- nsave*nskip+nburn
    
    p <- ns <- rep(0,L)
    v <- rep(1/L,L)
    v[L] <- 1
    prop <- prob <- matrix(0, nrow = n, ncol = L)
    
    z <- matrix(0, nrow = nsim, ncol = n)
    z[1,] <- rep(1,n)
    
    P <- Mu <- Mu1 <- Sigma2 <- Sigma21 <- matrix(0, nrow = nsim, ncol = L)
    Mu[1,] <- rep(mean(yt), L)
    Sigma2[1,] <- rep(var(yt), L)
    
    alpha <- numeric(nsim)
    alpha[1] <- 1
    
    for(i in 2:nsim) {
        cumv <- cumprod(1-v)
        p[1] <- v[1]
        for(l in 2:L){
            p[l] <- v[l]*cumv[l-1]
        }
        
        for(l in 1:L){
            prop[,l] <- p[l]*dnorm(yt, mean = Mu[i-1,l], sd = sqrt(Sigma2[i-1,l]))
        }
        prob <- prop/apply(prop,1,sum)
        for(j in 1:n){
            z[i,j] <- sample(1:L, size = 1, prob = prob[j, ])
        }
        P[i,] <- p
        
        for(l in 1:L){
            ns[l] <- length(which(z[i,] == l))
        }
        
        for(l in 1:(L-1)){
            v[l] <- rbeta(1, 1+ns[l], alpha[i-1] + sum(ns[(l+1):L]))
        }
        
        alpha[i] <- rgamma(1, shape = aalpha + L, balpha - sum(log(v[1:(L-1)])))
        
        for(l in 1:L){
            varmu <- 1/((1/S0) + (ns[l]/Sigma2[i-1,l]))
            meanmu <- ((sum(yt[z[i,] == l])/Sigma2[i-1,l]) + (m0/S0))/((1/S0) + (ns[l]/Sigma2[i-1,l]))
            Mu1[i,l] <- Mu[i,l] <- rnorm(1, mean = meanmu, sd = sqrt(varmu))
            if(standardise){
                Mu1[i,l] <- sd(y)*Mu[i,l]+mean(y)
            }
            
            Sigma21[i,l] <- Sigma2[i,l] <- 1/rgamma(1, a+ns[l]/2, b+0.5*sum((yt[z[i,]==l]-Mu[i,l])^2))
            if(standardise){
                Sigma21[i,l] <- var(y)*Sigma2[i,l]
            }
        }
    }
    
    res <- list()
    res$z <- z[seq(nburn+1, nsim, by = nskip),]
    res$P <- P[seq(nburn+1, nsim, by = nskip),]
    res$Mu <- Mu1[seq(nburn+1, nsim, by = nskip),]
    res$Sigma2 <- Sigma21[seq(nburn+1, nsim, by = nskip),]
    return(res)
}
