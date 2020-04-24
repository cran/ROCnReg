pooledROC.BB <-
function(marker, group, tag.h, data, p = seq(0, 1, l = 101), B = 5000, pauc = pauccontrol()) {
    
    pauc <- do.call("pauccontrol", pauc)
    
    # Obtain the marker in healthy and diseased
    yh <- data[,marker][data[,group] == tag.h]
    yd <- data[,marker][data[,group] != tag.h]
    
    # Missing data
    omit.h <- is.na(yh)
    omit.d <- is.na(yd)

    yh.wom <- yh[!omit.h]
    yd.wom <- yd[!omit.d]

    n1 <- length(yd.wom)
    n0 <- length(yh.wom)
    
    np <- length(p)
    
    u <- matrix(0, nrow = n1, ncol = B)
    u1 <- matrix(0, nrow = n0, ncol = B)

    weights.h <- matrix(0, nrow = n0, ncol = B)
    weights.d <- matrix(0, nrow = n1, ncol = B)
    
    rocbbpool <- matrix(0, nrow = np, ncol = B)
    aucbbpool <- numeric(B)
    
    if(pauc$compute) {
        paucbbpool <- numeric(B)
    }
    
    for(l in 1:B) {
        q <- rexp(n0, 1)
        weights.h[,l] <- q/sum(q)
        
        q1 <- rexp(n1, 1)
        weights.d[,l] <- q1/sum(q1)
        
        #u[,l] <- apply(outer(yh.wom, yd.wom, ">"), 2, weighted.mean, w = weights.h[,l])
        #rocbbpool[,l] <- apply(outer(u[,l], p, "<="), 2, weighted.mean, w = weights.d[,l])
        #I was unable to find the ewcdf function
        u[,l] <- 1 - ewcdf(yh.wom, weights.h[,l])(yd.wom)
        rocbbpool[,l] <- ewcdf(u[,l], weights.d[,l])(p)
        aucbbpool[l] <- sum((1-u[,l])* weights.d[,l])

		if(pauc$compute) {
			if(pauc$focus == "FPF") {
				paucbbpool[l] <- sum((pauc$value - pmin(u[,l], pauc$value))*weights.d[,l])
			} else {
				u1[,l] <- 1 - ewcdf(yd.wom, weights.d[,l])(yh.wom)
				paucbbpool[l] <- sum((pmax(u1[,l], pauc$value)-pauc$value)*weights.h[,l])
			}
		} 
        
    }
    poolROC <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
    poolROC[,1] <- apply(rocbbpool, 1, mean)
    poolROC[,2] <- apply(rocbbpool, 1, quantile, 0.025)
    poolROC[,3] <- apply(rocbbpool, 1, quantile, 0.975)
    
    res <- list()
    res$call <- match.call()
    res$marker <- list(h = yh, d = yd)
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$p <- p
    res$ROC <- poolROC
    AUC <- c(mean(aucbbpool), quantile(aucbbpool,c(0.025,0.975)))
    names(AUC) <- c("est","ql", "qh")
    res$AUC <- AUC
    res$weights <- list(h = weights.h, d = weights.d)
    if(pauc$compute) {
        pAUC <- c(mean(paucbbpool), quantile(paucbbpool, c(0.025,0.975)))
        names(pAUC) <- c("est","ql", "qh")
        res$pAUC <- if(pauc$focus == "FPF") {
        	pAUC/pauc$value
        } else {
        	pAUC/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }	
    class(res) <- c("pooledROC.BB", "pooledROC")
    return(res)
}
