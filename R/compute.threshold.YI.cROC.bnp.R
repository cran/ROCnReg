compute.threshold.YI.cROC.bnp <-
function(object, newdata) {
    if(class(object)[1] != "cROC.bnp") {
        stop(paste0("This function can not be used for this object class: ", class(object)[1]))
    }
    names.cov.h <- all.vars(object$fit$h$formula)[-1]
    names.cov.d <- all.vars(object$fit$d$formula)[-1]
    names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])
    
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    
    if(missing(newdata)) {
        newdata <- cROCData(object$data, names.cov, object$group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop=FALSE])
    }
    
    # Compute F_D|X and F_{\bar{D}}|X
    X0p <- predict(object$fit$h$mm, newdata = newdata)$X
    X1p <- predict(object$fit$d$mm, newdata = newdata)$X
    
    Lh <- object$prior$h$L
    Ld <- object$prior$d$L
    
    if(Lh > 1){
        Beta0 <- object$fit$h$beta
        Sigma0 <- object$fit$h$sd
        P0 <- object$fit$h$probs
    }
    
    if(Lh == 1){
        Beta0 <- object$fit$h$beta
        Sigma0 <- object$fit$h$sd
        P0 <- NULL
    }
    
    if(Ld > 1){
        Beta1 <- object$fit$d$beta
        Sigma1 <- object$fit$d$sd
        P1 <- object$fit$d$probs
    }
    
    if(Ld == 1){
        Beta1 <- object$fit$d$beta
        Sigma1 <- object$fit$d$sd
        P1 <- NULL
    }
    
    nsimf <- nrow(Beta0)
    npred <- nrow(newdata)
    
    y0 <- object$data_model$y$h
    y1 <- object$data_model$y$d
    n0 <- length(y0)
    n1 <- length(y1)
    
    grid  <- seq(min(c(y0, y1), na.rm = TRUE) - 1, max(c(y0, y1), na.rm = TRUE) + 1, length = max(500, c(n0,n1)))
    
    ngrid <- length(grid)
    
    F0 <- F1 <- array(0, c(ngrid, nsimf, npred))
    thresholds.s <- YI.s <- TPF.s <- FPF.s <- matrix(0, ncol = nsimf, nrow = npred)
    
    for(k in 1:nsimf) {
        if(Ld == 1 & Lh == 1){
            mu.h <- X0p%*%Beta0[k,]
            mu.d <- X1p%*%Beta1[k,]
            
            for(l in 1:npred){
                F0[,k,l] <- pnorm(grid, mu.h[l], Sigma0[k])
                F1[,k,l] <- pnorm(grid, mu.d[l], Sigma1[k])
                
                dif <- abs(F0[,k,l] -  F1[,k,l])
                thresholds.s[l,k] <- mean(grid[which(dif == max(dif))])
                YI.s[l,k] <- max(dif)
                
                TPF.s[l,k] <- 1 - pnorm(thresholds.s[l,k], mu.d[l], Sigma1[k])
                FPF.s[l,k] <- 1 - pnorm(thresholds.s[l,k], mu.h[l], Sigma0[k])
            }
        }
        
        if(Ld == 1 & Lh > 1){
            mu.d <- X1p%*%Beta1[k,]
            mu.h <- X0p%*%t(Beta0[k,,])
            for(l in 1:npred){
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = Sigma0[k,], w = P0[k,])
                F0[,k,l] <- pnorMix(grid, aux0)
                F1[,k,l] <- pnorm(grid, mu.d[l], Sigma1[k])
                
                dif <- abs(F0[,k,l] -  F1[,k,l])
                thresholds.s[l,k] <- mean(grid[which(dif == max(dif))])
                YI.s[l,k] <- max(dif)
                
                TPF.s[l,k] <- 1 - pnorm(thresholds.s[l,k], mu.d[l], Sigma1[k])
                FPF.s[l,k] <- 1 - pnorMix(thresholds.s[l,k], aux0)
            }
        }
        
        if(Ld > 1 & Lh == 1){
            mu.h <- X0p%*%Beta0[k,]
            mu.d <- X1p%*%t(Beta1[k,,])
            for(l in 1:npred){
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = Sigma1[k,], w = P1[k,])
                F0[,k,l] <- pnorm(grid, mu.h[l], Sigma0[k])
                F1[,k,l] <- pnorMix(grid, aux1)
                
                dif <- abs(F0[,k,l] -  F1[,k,l])
                thresholds.s[l,k] <- mean(grid[which(dif == max(dif))])
                YI.s[l,k] <- max(dif)
                
                TPF.s[l,k] <- 1 - pnorMix(thresholds.s[l,k], aux1)
                FPF.s[l,k] <- 1 - pnorm(thresholds.s[l,k], mu.h[l], Sigma0[k])
            }
        }
        
        if(Ld > 1 & Lh > 1){
            mu.h <- X0p%*%t(Beta0[k,,])
            mu.d <- X1p%*%t(Beta1[k,,])
            for(l in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = Sigma0[k,], w = P0[k,])
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = Sigma1[k,], w = P1[k,])
                
                F0[,k,l] <- pnorMix(grid, aux0)
                F1[,k,l] <- pnorMix(grid, aux1)
                
                difbb <- abs(F0[,k,l] -  F1[,k,l])
                thresholds.s[l,k] <- mean(grid[which(difbb == max(difbb))])
                YI.s[l,k] <- max(difbb)
                
                TPF.s[l,k] <- 1 - pnorMix(thresholds.s[l,k], aux1)
                FPF.s[l,k] <- 1 - pnorMix(thresholds.s[l,k], aux0)
            }
        }
    }
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$thresholds <- cbind(est = apply(thresholds.s, 1, mean), ql = apply(thresholds.s, 1, quantile, 0.025), qh = apply(thresholds.s, 1, quantile, 0.975))
    res$YI <- cbind(est = apply(YI.s, 1, mean), ql = apply(YI.s, 1, quantile, 0.025), qh = apply(YI.s, 1, quantile, 0.975))
    res$FPF <- cbind(est = apply(FPF.s, 1, mean), ql = apply(FPF.s, 1, quantile, 0.025), qh = apply(FPF.s, 1, quantile, 0.975))
    res$TPF <- cbind(est = apply(TPF.s, 1, mean), ql = apply(TPF.s, 1, quantile, 0.025), qh = apply(TPF.s, 1, quantile, 0.975))
    res
}
