compute.threshold.YI.pooledROC.dpm <-
function(object) {
    if(class(object)[2] != "pooledROC.dpm") {
        stop(paste0("This function can not be used for this object class: ", class(object)[2]))
    }
    
    y0 <- object$marker$h[!object$missing.ind$h]
    y1 <- object$marker$d[!object$missing.ind$d]
    
    p0 <- object$fit$h$P
    mu0 <- object$fit$h$Mu
    sigma02 <- object$fit$h$Sigma2
    p1 <- object$fit$d$P
    mu1 <- object$fit$d$Mu
    sigma12 <- object$fit$d$Sigma2
    
    if(is.null(p0) & is.null(p1)) {
        niter <- length(mu0)
    } else if(is.null(p0) & !is.null(p1)) {
        niter <- nrow(p1)
    } else if (!is.null(p0) & is.null(p1)) {
        niter <- nrow(p0)
    } else {
        niter <- nrow(p1)
    }
    
    grid <- seq(min(c(y0,y1))-1, max(c(y0,y1))+1, len = max(c(length(y0), length(y1)), 500))
    ngrid <- length(grid)
    
    #grid <- sort(unique(c(y0, y1)))
    #ngrid <- length(grid)
    
    thresholds.s <- YI.s <- FPF.s <- TPF.s <- numeric(niter)
    
    for(l in 1:niter){
        if(is.null(p0) & is.null(p1)){
            F0 <- pnorm(grid, mean = mu0[l], sd = sqrt(sigma02[l]))
            F1 <- pnorm(grid, mean = mu1[l], sd = sqrt(sigma12[l]))
            
            dif <- F0 - F1
            
            thresholds.s[l] <- mean(grid[which(dif == max(dif))])
            YI.s[l] <- max(dif)
            
            TPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = mu1[l], sd = sqrt(sigma12[l]))
            FPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = mu0[l], sd = sqrt(sigma02[l]))
        } else if(is.null(p0) & !is.null(p1)){
            aux1 <- norMix(mu = mu1[l,], sigma = sqrt(sigma12[l,]), w = p1[l,])
            
            F0 <- pnorm(grid, mean = mu0[l], sd = sqrt(sigma02[l]))
            F1 <- pnorMix(grid, aux1)
            
            dif <- F0 - F1
            
            thresholds.s[l] <- mean(grid[which(dif == max(dif))])
            YI.s[l] <- max(dif)
            
            TPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux1)
            FPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = mu0[l], sd = sqrt(sigma02[l]))
            
        } else if (!is.null(p0) & is.null(p1)){
            aux0 <- norMix(mu = mu0[l,], sigma = sqrt(sigma02[l,]), w = p0[l,])
            
            F1 <- pnorm(grid, mean = mu1[l], sd = sqrt(sigma12[l]))
            F0 <- pnorMix(grid, aux0)
            
            dif <- F0 - F1
            
            thresholds.s[l] <- mean(grid[which(dif == max(dif))])
            YI.s[l] <- max(dif)
            
            TPF.s[l] <- 1 - pnorm(thresholds.s[l], mean = mu1[l], sd = sqrt(sigma12[l]))
            FPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux0)
        } else {
            aux0 <- norMix(mu = mu0[l,], sigma = sqrt(sigma02[l,]), w = p0[l,])
            aux1 <- norMix(mu = mu1[l,], sigma = sqrt(sigma12[l,]), w = p1[l,])
            
            F0 <- pnorMix(grid, aux0)
            F1 <- pnorMix(grid, aux1)
            
            dif <- F0 - F1
            
            thresholds.s[l] <- mean(grid[which(dif == max(dif))])
            YI.s[l] <- max(dif)
            TPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux1)
            FPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux0)
        }
        
    }
    
    thresholds <- c(mean(thresholds.s), quantile(thresholds.s, c(0.025,0.975)))
    YI <- c(mean(YI.s), quantile(YI.s, c(0.025,0.975)))
    FPF <- c(mean(FPF.s), quantile(FPF.s, c(0.025,0.975)))
    TPF <- c(mean(TPF.s), quantile(TPF.s, c(0.025,0.975)))
    
    names(thresholds) <- names(YI) <- names(FPF) <- names(TPF) <- c("est","ql", "qh")
    
    res <- list()
    res$call <- match.call()
    res$thresholds <- thresholds
    res$YI <- YI
    res$FPF <- FPF
    res$TPF <- TPF
    res
}
