inf_criteria <-
function(y, X, res){
    n <- length(y)
    
    if(is.null(res$P)){
        L <- 1
    } else{
        L <- ncol(res$Beta)
    }
    
    if(L > 1){
        p <- res$P
    }
    if(L == 1){
        p <- NULL
    }
    
    beta <- res$Beta
    sigma2 <- res$Sigma2
    niter <- nrow(beta)
    
    if(L == 1){
        term <- matrix(0, nrow = niter, ncol = n)
        for(k in 1:niter) {
            term[k,] <- dnorm(y, mean = X%*%beta[k,], sd = sqrt(sigma2[k]))
        }        
    }
    
    if(L > 1){
        term_1 <- array(0, c(niter, L, n))
        term <- matrix(0, nrow = niter, ncol = n)
        
        for(i in 1:n) {
            for(l in 1:L) {
                term_1[,l,i] <- p[,l]*dnorm(y[i], mean = c(X[i,]%*%t(beta[,l,])), sd = sqrt(sigma2[,l]))
            }
            term[,i] <- apply(term_1[,,i], 1, function(x) sum(x))
        }
    }
    
    term
}
