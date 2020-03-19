compute.threshold.FPF.AROC.bnp <-
function(object, newdata, FPF = 0.5) {
    if(class(object)[2] != "AROC.bnp") {
        stop(paste0("This function cannot be used for this object class: ", class(object)[2]))
    }
    
    # Newdata
    names.cov <- all.vars(object$fit$formula)[-1]
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
    stop("Newdata must be a data frame")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
    stop("Not all needed variables are supplied in newdata")
    
    if(missing(newdata)) {
        newdata <- cROCData(object$data, names.cov, object$group)
    } else {
        newdata <- na.omit(newdata[,names.cov, drop = FALSE])
    }
    
    thresholds <- compute.threshold.FPF.bnp(object = object$fit, newdata = newdata, FPF = FPF)
    np <- dim(thresholds)[1]
    ncov <- dim(thresholds)[2]
    
    thresholdsm <- thresholdsl <- thresholdsh <- matrix(0, nrow = np, ncol = ncov)
    rownames(thresholdsm) <- rownames(thresholdsl) <- rownames(thresholdsh) <- FPF
    
    thresholdsm <- apply(thresholds, c(1,2), mean)
    thresholdsl <- apply(thresholds, c(1,2), quantile, 0.025)
    thresholdsh <- apply(thresholds, c(1,2), quantile, 0.975)
    
    # Organised results as desired
    thresholds.ret <- vector("list", length(FPF))
    names(thresholds.ret) <- FPF
    for(i in 1:length(FPF)){
        thresholds.ret[[i]] <- matrix(ncol = 3, nrow = ncov)
        thresholds.ret[[i]][,1] <- thresholdsm[i,] 
        thresholds.ret[[i]][,2] <- thresholdsl[i,]
        thresholds.ret[[i]][,3] <- thresholdsh[i,]
        colnames(thresholds.ret[[i]]) <- c("est", "ql", "qh")
    }


    res <- list()
    res$newdata <- newdata
    #res$thresholds.est <- thresholdsm
    #res$thresholds.ql <- thresholdsl
    #res$thresholds.qh <- thresholdsh
    res$thresholds <- thresholds.ret
    res$FPF <- FPF
    res
}
