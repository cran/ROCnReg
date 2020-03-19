bddp <-
function(y, X, prior, mcmc, standardise = TRUE) {
    yt <- y
    if(standardise == TRUE) {
        yt <- (y-mean(y))/sd(y)
        #yt <- y/sd(y)
    }
    n <- length(y)
    k <- ncol(X)
    
    m <- prior$m0
    S <- prior$S0
    nu <- prior$nu
    psi <- prior$Psi
    a <- prior$a
    b <- prior$b
    aalpha <- prior$aalpha
    balpha <- prior$balpha
    L <- prior$L
    
    nburn <- mcmc$nburn
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nsim <- nburn + nsave*nskip
    
    p <- ns <- rep(0,L)
    v <- rep(1/L,L)
    v[L] <- 1
    
    z <- matrix(0, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
    z[1,] <- rep(1,n)
    
    beta <- matrix(0, nrow = L, ncol = k)
    aux <- try(solve(t(X)%*%X)%*%t(X)%*%yt, silent = TRUE)
    if(!inherits(aux, "try-error")) {
        for(l in 1:L) {
            beta[l,] <- aux
        }
    }
    
    tau <- rep(1/var(yt),L)
    prop <- prob <- matrix(0, nrow = n, ncol = L)
    
    P <- Tau <- Sigma2 <- matrix(0, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
    Beta <- Beta1 <- array(0,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
    Beta[1,,] <- beta
    Tau[1,] <- tau
    
    mu <- matrix(0, nrow = nsim, ncol = k)
    Sigmainv <- array(0, c(nsim,k,k))
    mu[1,] <- mvrnorm(1, mu = m, Sigma = S) #rmvn(1, mu = m, sigma = S)
    Sigmainv[1,,] <- rWishart(1, df = nu, solve(nu*psi))
    
    alpha <- numeric(nsim)
    alpha[1]<-1
    
    for(i in 2:nsim) {
        cumv <- cumprod(1-v)
        p[1] <- v[1]
        for(l in 2:L) {
            p[l] <- v[l]*cumv[l-1]
        }
        for(l in 1:L) {
            prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,],sd=sqrt(1/tau[l]))
        }
        
        prob <- prop/apply(prop,1,sum)
        for(j in 1:n){
            z[i,j] <- sample(1:L, size = 1, prob = prob[j, ])
        }
        P[i,] <- p
        
        for(l in 1:L) {
            ns[l] <- length(which(z[i,]==l))
        }
        
        for(l in 1:(L-1)) {
            v[l] <- rbeta(1,1+ns[l],alpha[i-1]+sum(ns[(l+1):L]))
        }
        
        alpha[i] <- rgamma(1,shape=aalpha+L,balpha-sum(log(v[1:(L-1)])))
        
        for(l in 1:L) {
            tX  <- matrix(t(X[z[i,]==l, ]), nrow = k, ncol = ns[l])
            V <- solve(Sigmainv[i-1,,]+tau[l]*tX%*%X[z[i,] ==l,])
            mu1 <- V%*%(Sigmainv[i-1,,]%*%mu[i-1,]+tau[l]*tX%*%yt[z[i,] ==l])
            Beta1[i,l,] <- Beta[i,l,] <- beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V) #rmvn(1,mu=mu1,sigma=V)
            if (standardise == TRUE) {
                Beta1[i,l,1] <- sd(y)*Beta[i,l,1] + mean(y)
                if(k > 1) {
                    Beta1[i,l,2:k] <- sd(y)*Beta[i,l,2:k]
                }
            }
            #Tau[i,l] <- tau[l] <- rgamma(1,shape = a + (ns[l]/2), rate = b + 0.5*(t(yt[z[i,]==l]-X[z[i,] ==l,]%*%beta[l,])%*%(yt[z[i,] ==l]-X[z[i,] ==l,]%*%beta[l,])))            
            Tau[i,l] <- tau[l] <- rgamma(1,shape = a + (ns[l]/2), rate = b + 0.5*(t(yt[z[i,]==l]-X[z[i,] ==l,]%*%t(beta[l,,drop=FALSE]))%*%(yt[z[i,] ==l]-X[z[i,] ==l,]%*%t(beta[l,,drop=FALSE]))))
            
            Sigma2[i,l] = 1/Tau[i,l]
            if (standardise == TRUE) {
                Sigma2[i,l] <- var(y)*(1/Tau[i,l])
            }
        }
        
        Vaux <- solve(solve(S)+L*Sigmainv[i-1,,])
        if(k == 1) {
            meanmu <- Vaux%*%(solve(S)%*%m+Sigmainv[i-1,,]%*%sum(Beta[i,,]))
        } else {
            meanmu <- Vaux%*%(solve(S)%*%m+Sigmainv[i-1,,]%*%t(t(apply(Beta[i,,],2,sum))))
        }
        mu[i,] <- mvrnorm(1, mu = meanmu, Sigma = Vaux) #rmvn(1,mu = meanmu, sigma = Vaux)
        
        Vaux1 <- 0
        for(l in 1:L) {
            Vaux1 <- Vaux1+(Beta[i,l,]-mu[i,])%*%t((Beta[i,l,]-mu[i,]))
        }
        Sigmainv[i,,] <- rWishart(1,nu+L,solve(nu*psi+Vaux1))
    }
    res <- list()
    res$z <- z[seq(nburn+1, nsim, by = nskip),]
    res$P <- P[seq(nburn+1, nsim, by = nskip),]
    res$Beta <- Beta1[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
    res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
    res
}
