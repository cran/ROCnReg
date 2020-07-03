dpm <-
function(y, prior, mcmc, standardise = FALSE) {
    
    multinom <- function(probs) {
        probs <- t(apply(probs,1,cumsum)) 
        res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
        return(res)  
    }

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

    prop <- matrix(NA_real_, nrow = n, ncol = L)

    z <- matrix(NA_real_, nrow = nsim, ncol = n)
    z_tmp <- vector(length = n)

    z[1,] <- rep(1,n)
    
    #P <- Mu <- Mu1 <- Sigma2 <- Sigma21 <- matrix(0, nrow = nsim, ncol = L)
    P <- Mu <- Sigma2 <- matrix(NA_real_, nrow = nsim, ncol = L)

    #Mu_tmp <- Mu1_tmp <- Sigma2_tmp <- Sigma21_tmp <- vector(length = L)
    Mu_tmp <- Sigma2_tmp  <- vector(length = L)

    Mu[1,] <- rep(mean(yt), L)
    Sigma2[1,] <- rep(var(yt), L)
    
    Mu_tmp <- Mu[1,]
    Sigma2_tmp <- Sigma2[1,]

    alpha <- numeric(nsim)
    alpha[1] <- 1
    
    for(i in 2:nsim) {
        Sigma2_tmp <- Sigma2[i-1,]

        cumv <- cumprod(1-v)
        p[1] <- v[1]
        p[2:L] <- v[2:L]*cumv[1:(L-1)]
        
        #for(l in 2:L){
        #    p[l] <- v[l]*cumv[l-1]
        #}
        
        for(l in 1:L){
            prop[,l] <- p[l]*dnorm(yt, mean = Mu_tmp[l], sd = sqrt(Sigma2_tmp[l]))
        }
        #prob <- prop/apply(prop,1,sum)
        prob <- prop/rowSums(prop)

        #for(j in 1:n){
        #    z[i,j] <- sample(1:L, size = 1, prob = prob[j, ])
        #}

        #z_tmp <- apply(prob, 1, function(x, L)  sample(1:L, 1, prob = x), L = L)
        #z_tmp <- apply(prob, 1, function(x, L)  sample.int(L, 1, prob = x), L = L)
        
        z_tmp <- multinom(prob)
        #for(l in 1:L){
        #    #ns[l] <- length(which(z[i,] == l))
        #    ns[l] <- sum(z[i,] == l)
        #}
        ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
        yt_z_l <- sapply(1:L, function(x, v, y) sum(y[v == x]), v = z_tmp, y = yt)

        #for(l in 1:(L-1)){
        #    v[l] <- rbeta(1, 1 + ns[l], alpha[i-1] + sum(ns[(l+1):L]))
        #}
        v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha[i-1] + rev(cumsum(rev(ns[-1]))))
        
        alpha[i] <- rgamma(1, shape = aalpha + L, balpha - sum(log(v[1:(L-1)])))

        varmu <- 1/((1/S0) + (ns/Sigma2_tmp))
        meanmu <- ((yt_z_l/Sigma2_tmp) + (m0/S0))/((1/S0) + (ns/Sigma2_tmp))
        Mu_tmp <- rnorm(L, mean = meanmu, sd = sqrt(varmu))
        #if(standardise){
        #    Mu1_tmp <- sd(y)*Mu_tmp + mean(y)
        #}
        yt_z_l_mu <- sapply(1:L, function(x, v, y, mu) sum((y[v == x] - mu[x])^2), v = z_tmp, y = yt, mu = Mu_tmp)

        Sigma2_tmp <- 1/rgamma(L, a + ns/2, b + 0.5*yt_z_l_mu)
        #if(standardise){
        #    Sigma21_tmp <- var(y)*Sigma2_tmp
        #}

        P[i,] <- p
        z[i,] <- z_tmp
        Mu[i,] <- Mu_tmp
        Sigma2[i,] <- Sigma2_tmp
        #Sigma21[i,] <- Sigma21_tmp

        #for(l in 1:L){
        #    varmu <- 1/((1/S0) + (ns[l]/Sigma2[i-1,l]))
        #    meanmu <- ((sum(yt[z[i,] == l])/Sigma2[i-1,l]) + (m0/S0))/((1/S0) + (ns[l]/Sigma2[i-1,l]))
        #    Mu1[i,l] <- Mu[i,l] <- rnorm(1, mean = meanmu, sd = sqrt(varmu))
        #    if(standardise){
        #        Mu1[i,l] <- sd(y)*Mu[i,l]+mean(y)
        #    }
            
        #    Sigma21[i,l] <- Sigma2[i,l] <- 1/rgamma(1, a+ns[l]/2, b+0.5*sum((yt[z[i,]==l]-Mu[i,l])^2))
        #    if(standardise){
        #        Sigma21[i,l] <- var(y)*Sigma2[i,l]
        #    }
        #}
    }
    if(standardise){
        Mu <- sd(y)*Mu + mean(y)
        Sigma2 <- var(y)*Sigma2
    }

    res <- list()
    res$z <- z[seq(nburn+1, nsim, by = nskip),]
    res$P <- P[seq(nburn+1, nsim, by = nskip),]
    res$Mu <- Mu[seq(nburn+1, nsim, by = nskip),]
    res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
    return(res)
}
