predictive.checks.cROC.bnp <-
function(object, statistics = c("min","max","kurtosis","skewness"), ndensity = 512, devnew = TRUE) {
    if(class(object)[1] != "cROC.bnp") {
        stop(paste0("This function can not be used for this object class: ", class(object)[1]))
    }
    
    # Nondiseased group
    y0 <- object$data_model$y$h
    n0 <- length(y0)
    nrep0 <- nrow(object$fit$h$beta)
    yrep0 <- matrix(0, nrow = n0, ncol = nrep0)
    
    if(object$prior$h$L > 1){
        aux <- t(apply(object$fit$h$probs[1:nrep0,], 1, function(x, n)  sample(1:length(x), n, replace = TRUE, prob = x), n = n0))
        for(l in 1:nrep0) {
            yrep0[,l] <- rnorm(n = n0, mean = colSums(t(object$data_model$X$h)*t(object$fit$h$beta[l,aux[l,],])), sd = object$fit$h$sd[l,aux[l,]])
        }
    }
    if(object$prior$h$L == 1){
        for(l in 1:nrep0) {
            yrep0[,l] <- rnorm(n = n0, mean = as.numeric(object$data_model$X$h%*%object$fit$h$beta[l,]), sd = object$fit$h$sd[l])
        }
    }
    predictive.checks.helper(y = y0, yrep = yrep0, object$marker, statistics = statistics, devnew = devnew, group = "H", ndensity = ndensity)
    
    # Diseased group
    y1 <- object$data_model$y$d
    n1 <- length(y1)
    nrep1 <- nrow(object$fit$d$beta)
    yrep1 <- matrix(0, nrow = n1, ncol = nrep1)
    
    if(object$prior$d$L > 1){
        aux <- t(apply(object$fit$d$probs[1:nrep1,], 1, function(x, n)  sample(1:length(x), n, replace = TRUE, prob = x), n = n1))
        for(l in 1:nrep1) {
            yrep1[,l] <- rnorm(n = n1, mean = colSums(t(object$data_model$X$d)*t(object$fit$d$beta[l,aux[l,],])), sd = object$fit$d$sd[l,aux[l,]])
        }
    }
    if(object$prior$d$L == 1){
        for(l in 1:nrep1) {
            yrep1[,l] <- rnorm(n = n1, mean = as.numeric(object$data_model$X$d%*%object$fit$d$beta[l,]), sd = object$fit$d$sd[l])
        }
    }
    if(devnew) dev.new()
    predictive.checks.helper(y = y1, yrep = yrep1, object$marker, statistics = statistics, devnew = devnew, group = "D", ndensity = ndensity)
    
    res <- list()
    res$yrep <- list(h = yrep0, d = yrep1)
    res$y <- list(h = y0, d = y1)
    invisible(res)
}
