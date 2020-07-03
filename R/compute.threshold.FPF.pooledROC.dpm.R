compute.threshold.FPF.pooledROC.dpm <-
function(object, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	doMCMCTH <- function(k, res0, res1, FPF) {
		p0 <- res0$P
		p1 <- res1$P

		if(is.null(p0) & is.null(p1)) {
            thresholds.s <- qnorm(1 - FPF, mean = res0$Mu[k], sd= sqrt(res0$Sigma2[k]))
            TPF.s <- 1 - pnorm(thresholds.s, mean = res1$Mu[k], sd = sqrt(res1$Sigma2[k]))
        } else if(is.null(p0) & !is.null(p1)){
            aux1 <- norMix(mu = res1$Mu[k,], sigma = sqrt(res1$Sigma2[k,]), w = p1[k,])
            thresholds.s <- qnorm(1 - FPF, mean = res0$Mu[k], sd= sqrt(res0$Sigma2[k]))
            TPF.s <- 1 - pnorMix(thresholds.s, aux1)
        } else if (!is.null(p0) & is.null(p1)){
            aux0 <- norMix(mu = res0$Mu[k,], sigma = sqrt(res0$Sigma2[k,]), w = p0[k,])
            thresholds.s <- qnorMix(1 - FPF, aux0)
            TPF.s <- 1 - pnorm(thresholds.s, mean = res1$Mu[k], sd = sqrt(res1$Sigma2[k]))
        } else {
            aux0 <- norMix(mu = res0$Mu[k,], sigma = sqrt(res0$Sigma2[k,]), w = p0[k,])
            aux1 <- norMix(mu = res1$Mu[k,], sigma = sqrt(res1$Sigma2[k,]), w = p1[k,])
            thresholds.s <- qnorMix(1 - FPF, aux0)
            TPF.s <- 1 - pnorMix(thresholds.s, aux1)
        }
        res <- list()
        res$thresholds.s <- thresholds.s
        res$TPF.s <- TPF.s
        res

	}

	if(class(object)[1] != "pooledROC.dpm") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}

	parallel <- match.arg(parallel)

	if(object$mcmc$nsave > 0) {
        do_mc <- do_snow <- FALSE
        if (parallel != "no" && ncpus > 1L) {
            if (parallel == "multicore") {
                do_mc <- .Platform$OS.type != "windows"
            } else if (parallel == "snow") {
                do_snow <- TRUE
            }
            if (!do_mc && !do_snow) {
                ncpus <- 1L
            }       
            loadNamespace("parallel") # get this out of the way before recording seed
        }
        # Seed
        #if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
        #seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

        # Apply function
        resBoot <- if (ncpus > 1L && (do_mc || do_snow)) {
                if (do_mc) {
                    parallel::mclapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, FPF = FPF, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, FPF = FPF)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, FPF = FPF)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, FPF = FPF)
            }

        resBoot <- simplify2array(resBoot)    
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        TPF.s <- simplify2array(resBoot["TPF.s",])
        if(length(FPF) == 1) {            	
        	thresholds.s <- matrix(thresholds.s, nrow = 1)
        	TPF.s <- matrix(TPF.s, nrow = 1)
        }

    } else {
        stop("nsave should be larger than zero.")
    }

    np <- length(FPF)
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
