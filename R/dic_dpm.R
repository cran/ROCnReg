dic_dpm <-
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
		logtermsum <- apply(logterm, 1, sum)

		deviance <- -2*mean(logtermsum)
		muhat <- mean(mu)
		sigma2hat <- mean(sigma2)
		d <- -2*sum( dnorm(y, muhat, sqrt(sigma2hat), log=TRUE))
		pd <- deviance - d
		DIC <- deviance + pd

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
		deviance <- -2 * mean(apply(logtermsum,1,sum))
		pd <- deviance + 2 * sum(log(apply(exp(logtermsum),2,mean)))
		DIC <- deviance + pd
	}

	res <- list()
	res$pD <- pd
	res$DIC <- DIC
	return(res)
}
