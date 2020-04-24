compute.threshold.FPF.cROC.bnp <-
function(object, newdata, FPF = 0.5) {
    if(class(object)[1] != "cROC.bnp") {
        stop(paste0("This function cannot be used for this object class: ", class(object)[1]))
    }
    # Newdata
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
    
    thresholds <- compute.threshold.FPF.bnp(object = object$fit$h, newdata = newdata, FPF = FPF)
    
    np <- dim(thresholds)[1]
    ncov <- dim(thresholds)[2]
    
    thresholdsm <- thresholdsl <- thresholdsh <- matrix(0, nrow = np, ncol = ncov)
    rownames(thresholdsm) <- rownames(thresholdsl) <- rownames(thresholdsh) <- FPF
    
    thresholdsm <- apply(thresholds, c(1,2), mean)
    thresholdsl <- apply(thresholds, c(1,2), quantile, 0.025)
    thresholdsh <- apply(thresholds, c(1,2), quantile, 0.975)
    
    # Compute associated TPF
    Xp <- predict.design.matrix.bnp(object$fit$d$mm, newdata)$X
    
    ncov <- nrow(Xp)
    nrep <- nrow(object$fit$d$beta)
    np <- length(FPF)
    
    TPF <- array(0,c(np,ncov,nrep))
    if(object$prior$d$L > 1){
        for(inrep in 1:nrep) {
            mu.d <- Xp%*%t(object$fit$d$beta[inrep,,])
            for(incov in 1:ncov) {
                aux <- norMix(mu = c(mu.d[incov,]), sigma = object$fit$d$sd[inrep,], w = object$fit$d$probs[inrep,])
                TPF[,incov,inrep] <- 1-pnorMix(thresholds[,incov,inrep], aux)
            }
        }
    }
    
    if(object$prior$d$L == 1){
        for(inrep in 1:nrep) {
            mu.d <- Xp%*%object$fit$d$beta[inrep,]
            for(incov in 1:ncov) {
                TPF[,incov,inrep] <- 1-pnorm(thresholds[,incov,inrep], mu.d[incov], object$fit$d$sd[inrep])
            }
        }
    }
    
    TPFm <- apply(TPF, c(1,2), mean)
    TPFl <- apply(TPF, c(1,2), quantile, 0.025)
    TPFh <- apply(TPF, c(1,2), quantile, 0.975)

    # Organised results as desired
    TPF.ret <- thresholds.ret <- vector("list", length(FPF))
    names(TPF.ret) <- names(thresholds.ret) <- FPF
    for(i in 1:length(FPF)){
        TPF.ret[[i]] <- thresholds.ret[[i]] <- matrix(ncol = 3, nrow = ncov)

        thresholds.ret[[i]][,1] <- thresholdsm[i,] 
        thresholds.ret[[i]][,2] <- thresholdsl[i,]
        thresholds.ret[[i]][,3] <- thresholdsh[i,]

        TPF.ret[[i]][,1] <- TPFm[i,] 
        TPF.ret[[i]][,2] <- TPFl[i,]
        TPF.ret[[i]][,3] <- TPFh[i,]

        colnames(TPF.ret[[i]]) <- colnames(thresholds.ret[[i]]) <- c("est", "ql", "qh")
    }
        
    res <- list()
    res$newdata <- newdata
    #res$thresholds.est <- thresholdsm
    #res$thresholds.ql <- thresholdsl
    #res$thresholds.qh <- thresholdsh
    #res$TPF.est <- TPFm
    #res$TPF.ql <- TPFl
    #res$TPF.qh <- TPFh
    res$thresholds <- thresholds.ret
    res$TPF <- TPF.ret
    res$FPF <- FPF
    res
}
