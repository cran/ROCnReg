waic_dpm <-
function(y, res, L, termsum = NULL) {  
	if(L == 1) {
		n <- length(y)		
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- length(mu)

		if(is.null(termsum)) {
			termsum <- inf_criteria_dpm(y = y, res = res, L = L)
		}

		logterm <- log(termsum)
		lpd <- sum(log(apply(exp(logterm),2,mean)))
		p2 <- sum(apply(logterm,2,var))
		waic <- -2*(lpd-p2)
	} else {
		n <- length(y)	
		p <- res$P
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- nrow(p)

		if(is.null(termsum)) {
			termsum <- inf_criteria_dpm(y = y, res = res, L = L)
		}

		logtermsum <- log(termsum)
		lpd <- sum(log(apply(exp(logtermsum),2,mean)))
		p2 <- sum(apply(logtermsum,2,var))
		waic <- -2*(lpd-p2)
	}

	res <- list()
	res$pW <- p2
	res$WAIC <- waic
	return(res)
}
