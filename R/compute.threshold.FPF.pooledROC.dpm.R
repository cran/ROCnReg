compute.threshold.FPF.pooledROC.dpm <-
function(object, FPF = 0.5) {
	if(class(object)[2] != "pooledROC.dpm") {
		stop(paste0("This function can not be used for this object class: ", class(object)[2]))
	}

	p0 <- object$fit$h$P
	mu0 <- object$fit$h$Mu
	sigma02 <- object$fit$h$Sigma2
	p1 <- object$fit$d$P
	mu1 <- object$fit$d$Mu
	sigma12 <- object$fit$d$Sigma2

	if(is.null(p0) & is.null(p1)) {
		niter <- length(mu0)
	} else if(is.null(p0) & !is.null(p1)) {
		niter <- nrow(p1)
	} else if (!is.null(p0) & is.null(p1)) {
		niter <- nrow(p0)
	} else {
		niter <- nrow(p1)
	}

	np <- length(FPF)
	thresholds.s <- TPF.s <- matrix(0, nrow = np, ncol = niter)

	for(k in 1:niter) {
		 if(is.null(p0) & is.null(p1)){
            thresholds.s[,k] <- qnorm(1 - FPF, mean = mu0[k], sd= sqrt(sigma02[k]))
            TPF.s[,k] <- 1 - pnorm(thresholds.s[,k], mean = mu1[k], sd = sqrt(sigma12[k]))
        } else if(is.null(p0) & !is.null(p1)){
            aux1 <- norMix(mu = mu1[k,], sigma = sqrt(sigma12[k,]), w = p1[k,])
            thresholds.s[,k] <- qnorm(1 - FPF, mean = mu0[k], sd= sqrt(sigma02[k]))
            TPF.s[,k] <- 1 - pnorMix(thresholds.s[,k], aux1)
        } else if (!is.null(p0) & is.null(p1)){
            aux0 <- norMix(mu = mu0[k,], sigma = sqrt(sigma02[k,]), w = p0[k,])
            thresholds.s[,k] <- qnorMix(1 - FPF, aux0)
            TPF.s[,k] <- 1 - pnorm(thresholds.s[,k], mean = mu1[k], sd = sqrt(sigma12[k]))
        } else {
            aux0 <- norMix(mu = mu0[k,], sigma = sqrt(sigma02[k,]), w = p0[k,])
            aux1 <- norMix(mu = mu1[k,], sigma = sqrt(sigma12[k,]), w = p1[k,])
            thresholds.s[,k] <- qnorMix(1 - FPF, aux0)
            TPF.s[,k] <- 1 - pnorMix(thresholds.s[,k], aux1)
        }
	}

	thresholds <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
	rownames(thresholds) <- FPF

	thresholds[,1] <- apply(thresholds.s, 1, mean)
	thresholds[,2] <- apply(thresholds.s, 1, quantile, prob = 0.025)
	thresholds[,3] <- apply(thresholds.s, 1, quantile, prob = 0.975)

	TPF <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
	rownames(TPF) <- FPF

	TPF[,1] <- apply(TPF.s, 1, mean)
	TPF[,2] <- apply(TPF.s, 1, quantile, prob = 0.025)
	TPF[,3] <- apply(TPF.s, 1, quantile, prob = 0.975)

	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res
}
