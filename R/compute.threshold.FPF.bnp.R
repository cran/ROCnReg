compute.threshold.FPF.bnp <-
function(object, newdata, FPF = 0.5) {
    Xp <- predict.design.matrix.bnp(object$mm, newdata)$X
    
    ncov <- nrow(Xp)
    nrep <- nrow(object$beta)
    np <- length(FPF)
    if(is.null(object$probs)){
        L = 1
    } else{
        L = ncol(object$beta)
    }
    
    thresholds <- array(0,c(np,ncov,nrep))
    
    if(L == 1){
        for(inrep in 1:nrep) {
            for(incov in 1:ncov) {
                thresholds[,incov,inrep] <- qnorm(1-FPF, mean = Xp[incov,]%*%object$beta[inrep,], sd = object$sd[inrep])
            }
        }
    }
    
    if(L > 1){
        for(inrep in 1:nrep) {
            mu.h <- Xp%*%t(object$beta[inrep,,])
            for(incov in 1:ncov) {
                aux <- norMix(mu = c(mu.h[incov,]), sigma = object$sd[inrep,], w = object$probs[inrep,])
                thresholds[,incov,inrep] <- qnorMix(1-FPF, aux)
            }
        }
    }
    
    thresholds
}
