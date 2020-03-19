lpml_dpm <-
function(y, res, L, termsum = NULL){
	if(L == 1) {
		n <- length(y)
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- length(mu)
		
		if(is.null(termsum)) {
			termsum <- inf_criteria_dpm(y = y, res = res, L = L)
		}
		
		term1 <- 1/termsum
		omegabari <- apply(term1,2,mean)
		omegabari_1 <- sqrt(niter) * omegabari
		omegatilde <- matrix(0, nrow = niter, ncol = n)

		for(i in 1:n) {
			omegatilde[,i] <- pmin(term1[,i], omegabari_1[i])
		}

		sum_omegatilde <- apply(omegatilde, 2 ,sum)
		sum_term_omegatilde <- apply(termsum*omegatilde, 2, sum)
		cpo <- sum_term_omegatilde/sum_omegatilde

		lpml <- sum(log(cpo))

	} else {
		n <- length(y)
		p <- res$P
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- nrow(p)
		
		if(is.null(termsum)) {
		    termsum <- inf_criteria_dpm(y = y, res = res, L = L)
		}
	    
	    aux <- 1/termsum
		omegabari <- apply(aux, 2, mean)
		omegabari_1 <- sqrt(niter) * omegabari
		omegatilde <- matrix(0, nrow = niter, ncol = n)

		for(i in 1:n) {
			omegatilde[,i] <- pmin(aux[,i], omegabari_1[i])  
		}

		sum_omegatilde <- apply(omegatilde,2,sum)
		sum_term_omegatilde <- apply(termsum*omegatilde, 2, sum)
		cpo <- sum_term_omegatilde/sum_omegatilde

		lpml <- sum(log(cpo))
	}
	res <- list()
	res$cpo <- cpo  
	res$lpml <- lpml
	return(res)
}
