compute.threshold.YI.AROC.bnp <-
function(object, newdata) {
    if(class(object)[2] != "AROC.bnp") {
        stop(paste0("This function cannot be used for this object class: ", class(object)[2]))
    }
    names.cov <- all.vars(object$fit$formula)[-1]
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    
    if(missing(newdata)) {
        newdata <- cROCData(object$data, names.cov, object$group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop=FALSE])
    }
    p <- seq(0, 1, length = 500)
    np <- length(p)
    npred <- nrow(newdata)
    
    # Compute AROC
    if(object$prior$L > 1){
        Beta0 <- object$fit$beta
        Sigma0 <- object$fit$sd
        P0 <- object$fit$probs
    }
    if(object$prior$L == 1){
        Beta0 <- object$fit$beta
        Sigma0 <- object$fit$sd
        P0 <- NULL
    }
    
    nsimf <- nrow(Beta0)
    n1 <- length(object$data_model$y$d)
    
    up <- matrix(0, nrow = n1, ncol = nsimf)
    weights <- matrix(0, nrow = n1, ncol = nsimf)
    AROC <- matrix(0, nrow = np, ncol = nsimf)
    
    if(object$prior$L > 1){
        for(l in 1:nsimf) {
            up[,l] <- 1 - apply(t(P0[l,]*t(pnorm(object$data_model$y$d, mean = object$data_model$X$d%*%t(Beta0[l,,]), sd = rep(Sigma0[l,], each = n1)))),1, sum)
            aux1 <- rexp(n1,1)
            weights[,l] <- aux1/sum(aux1)
        }
        for(j in 1:np) {
            AROC[j,] <- colSums(weights*(up <= p[j]))
        }
    }
    
    if(object$prior$L == 1){
        for(l in 1:nsimf) {
            up[,l] = 1 - pnorm(object$data_model$y$d, mean = object$data_model$X$d%*%Beta0[l,], sd = Sigma0[l])
            aux1 <- rexp(n1,1)
            weights[,l] <- aux1/sum(aux1)
        }
        for(j in 1:np) {
            AROC[j,] <- colSums(weights*(up<=p[j]))
        }
    }
    
    # Compute YI and associated threshold values
    X0p <- predict(object$fit$mm, newdata = newdata)$X
    thresholds <- matrix(0, ncol = nsimf, nrow = npred)
    YI.s <- FPF.s <- vector(length = nsimf)
    for(l in 1:nsimf) {
        dif <- AROC[,l] -  p
        FPF.s[l] <- mean(p[which(dif == max(dif))])
        YI.s[l] <- max(dif)
        if(object$prior$L > 1){
            mu.h <- X0p%*%t(Beta0[l,,])
            for(k in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[k,]), sigma = Sigma0[l,], w = P0[l,])
                thresholds[k,l] <- qnorMix(1 - FPF.s[l], aux0)
            }
        }
        if(object$prior$L == 1){
            mu.h <- X0p%*%Beta0[l,]
            for(k in 1:npred) {
                thresholds[k,l] <- qnorm(1 - FPF.s[l], mu.h[k], Sigma0[l])
            }
        }
    }
    
    YI <- c(mean(YI.s), quantile(YI.s, c(0.025,0.975)))
    FPF <- c(mean(FPF.s), quantile(FPF.s, c(0.025,0.975)))
    names(YI) <- names(FPF) <- c("est","ql", "qh")
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$thresholds <- cbind(est = apply(thresholds, 1, mean), ql = apply(thresholds, 1, quantile, 0.025), qh = apply(thresholds, 1, quantile, 0.975))
    res$YI <- YI
    res$FPF <- FPF
    res
}
