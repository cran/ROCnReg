predictive.checks.AROC.bnp <-
function(object, statistics = c("min","max","kurtosis","skewness"), ndensity = 512, devnew = TRUE) {
    if(class(object)[1] != "AROC.bnp") {
        stop(paste0("This function can not be used for this object class: ", class(object)[1]))
    }
    
    y0 <- object$data_model$y$h
    n0 <- length(y0)
    nrep <- nrow(object$fit$beta)
    yrep <- matrix(0, nrow = n0, ncol = nrep)
    
    if(object$prior$L > 1){
        aux <- t(apply(object$fit$probs[1:nrep,], 1, function(x, n)  sample(1:length(x), n, replace = TRUE, prob = x), n = n0))
        for(l in 1:nrep) {
            yrep[,l] <- rnorm(n = n0, mean = colSums(t(object$data_model$X$h)*t(object$fit$beta[l,aux[l,],])), sd = object$fit$sd[l,aux[l,]])
        }
    }
    if(object$prior$L == 1){
        for(l in 1:nrep) {
            yrep[,l] <- rnorm(n = n0, mean = as.numeric(object$data_model$X$h%*%object$fit$beta[l,]), sd = object$fit$sd[l])
        }
    }
    
    predictive.checks.helper(y = y0, yrep = yrep, object$marker, statistics = statistics, devnew = devnew, group = "H", ndensity = ndensity)
    
    res <- list()
    res$yrep <- list(h = yrep)
    res$y <- list (h = y0)
    invisible(res)
}
