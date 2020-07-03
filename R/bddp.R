bddp <-
function(y, X, prior, mcmc, standardise = TRUE) {
    multinom <- function(prob) {
        probs <- t(apply(prob,1,cumsum)) 
        res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
        return(res)  
    }

    yt <- y
    if(standardise == TRUE) {
        yt <- (y-mean(y))/sd(y)
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
    
    p <- ns <- rep(0, L)
    v <- rep(1/L,L)
    v[L] <- 1
    
    z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
    z_tmp <- vector(length = n)

    z[1,] <- rep(1,n)
    
    beta <- matrix(NA_real_, nrow = L, ncol = k)
    aux <- ols.function(X, yt)$coeff
    if(!inherits(aux, "try-error")) {
        for(l in 1:L) {
            beta[l,] <- aux
        }
    }
    
    tau <- rep(1/var(yt),L)
    prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
    
    P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
    Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
    
    Beta[1,,] <- beta
    Tau[1,] <- tau
    
    #mu <- matrix(NA_real_, nrow = nsim, ncol = k)
    #Sigmainv <- array(NA_real_, c(nsim, k, k))
    
    #mu[1,] <- mvrnorm(1, mu = m, Sigma = S) #rmvn(1, mu = m, sigma = S)
    #Sigmainv[1,,] <- rWishart(1, df = nu, solve(nu*psi))
    
    mu <- mvrnorm(1, mu = m, Sigma = S)
    Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 

    alpha <- numeric(nsim)
    alpha[1] <- 1
    
    for(i in 2:nsim) {
        cumv <- cumprod(1-v)
        p[1] <- v[1]
        #for(l in 2:L) {
        #    p[l] <- v[l]*cumv[l-1]
        #}
        cumv <- cumprod(1-v)
        p[1] <- v[1]
        p[2:L] <- v[2:L]*cumv[1:(L-1)]

        for(l in 1:L) {
            prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
        }
        
        #prob <- prop/apply(prop,1,sum)
        prob <- prop/rowSums(prop)

        #for(j in 1:n){
        #    z[i,j] <- sample(1:L, size = 1, prob = prob[j, ])
        #}
        z_tmp <- multinom(prob)

        
        #for(l in 1:L) {
        #    ns[l] <- length(which(z_tmp == l))
        #}

        ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)

        #for(l in 1:(L-1)) {
        #    v[l] <- rbeta(1, 1 + ns[l], alpha[i-1] + sum(ns[(l+1):L]))
        #}
        v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
        
        alpha[i] <- rgamma(1, shape = aalpha + L, balpha - sum(log(v[1:(L-1)])))
        #Sigmainv_mu <- Sigmainv[i-1,,]%*%mu[i-1,]
        Sigmainv_mu <- Sigmainv%*%mu
        
        for(l in 1:L) {
            X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
            yt_l <- yt[z_tmp == l]
            #V <- solve(Sigmainv[i-1,,] + tau[l]*crossprod(X_l))
            V <- solve(Sigmainv + tau[l]*crossprod(X_l))
            mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
            beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)

            aux <- yt_l - X_l%*%beta[l,]
            tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
        }

        S_inv <- solve(S)
        #Vaux <- solve(S_inv + L*Sigmainv[i-1,,])
        Vaux <- solve(S_inv + L*Sigmainv)
        if(k == 1) {
            #meanmu <- Vaux%*%(S_inv%*%m + Sigmainv[i-1,,]%*%sum(beta))
            meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
        } else {
            #meanmu <- Vaux%*%(S_inv%*%m + Sigmainv[i-1,,]%*%colSums(beta))
            meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
        }
        #mu[i,] <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
        mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
        
        Vaux1 <- 0
        for(l in 1:L) {
            #Vaux1 <- Vaux1 + tcrossprod(beta[l,]-mu[i,])
            Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
        }
        
        #Sigmainv[i,,] <- rWishart(1, nu + L, solve(nu*psi + Vaux1))
        Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]

        P[i,] <- p
        z[i,] <- z_tmp
        Beta[i,,] <- beta
        Tau[i,] <- tau
    }

    if (standardise == TRUE) {
        Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
        if(k > 1) {
            Beta[,,2:k] <- sd(y)*Beta[,,2:k]
        }
        Sigma2 <- var(y)*(1/Tau)
    } else {
        Sigma2 <- 1/Tau
    }

    res <- list()
    res$z <- z[seq(nburn+1, nsim, by = nskip),]
    res$P <- P[seq(nburn+1, nsim, by = nskip),]
    res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
    res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
    res
}
