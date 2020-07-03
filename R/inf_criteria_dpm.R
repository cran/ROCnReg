inf_criteria_dpm <-
function(y, res, L) {
	if(L == 1) {
		n <- length(y)  
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- length(mu)

		termsum <- matrix(0, nrow = niter, ncol = n)
		#for(k in 1:niter) {
		#	termsum[k,] <- dnorm(y, mean = mu[k], sd = sqrt(sigma2[k]))
		#}
		for(i in 1:n) {
			termsum[,i] <- dnorm(y[i], mean = mu, sd = sqrt(sigma2))
		}
	} else {
		n <- length(y)
		p <- res$P
		mu <- res$Mu
		sigma2 <- res$Sigma2
		niter <- nrow(p)
		term <- array(0, c(niter, L, n))
	  	termsum <- matrix(0, nrow = niter, ncol = n)
	  	for(i in 1:n) {
	  		for(l in 1:L) {
	  			term[,l,i] <- p[,l]*dnorm(y[i], mean = mu[,l], sd = sqrt(sigma2[,l]))
	  		}
	  		termsum[,i] <- apply(term[,,i], 1, function(x) sum(x))
	  	}
	}
	termsum	
}
