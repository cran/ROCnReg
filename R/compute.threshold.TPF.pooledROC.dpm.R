compute.threshold.TPF.pooledROC.dpm <-
function(object, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	doMCMCTH <- function(k, res0, res1, TPF) {
		p0 <- res0$probs
		p1 <- res1$probs

		if(is.null(p0) & is.null(p1)) {
            thresholds.s <- qnorm(1 - TPF, mean = res1$mu[k], sd= res1$sd[k])
            FPF.s <- 1 - pnorm(thresholds.s, mean = res0$mu[k], sd = res0$sd[k])
        } else if(is.null(p0) & !is.null(p1)){
            aux0 <- norMix(mu = res0$mu[k,], sigma = res0$sd[k,], w = p0[k,])
            thresholds.s <- qnorm(1 - TPF, mean = res1$mu[k], sd= res1$sd[k])
            FPF.s <- 1 - pnorMix(thresholds.s, aux0)
        } else if (!is.null(p0) & is.null(p1)){
            aux1 <- norMix(mu = res1$mu[k,], sigma = res1$sd[k,], w = p1[k,])
            thresholds.s <- qnorMix(1 - TPF, aux1)
            FPF.s <- 1 - pnorm(thresholds.s, mean = res0$mu[k], sd = res0$sd[k])
        } else {
            aux1 <- norMix(mu = res1$mu[k,], sigma = res1$sd[k,], w = p1[k,])
            aux0 <- norMix(mu = res0$mu[k,], sigma = res0$sd[k,], w = p0[k,])
            thresholds.s <- qnorMix(1 - TPF, aux1)
            FPF.s <- 1 - pnorMix(thresholds.s, aux0)
        }
        res <- list()
        res$thresholds.s <- thresholds.s
        res$FPF.s <- FPF.s
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
                    parallel::mclapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, TPF = TPF, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, TPF = TPF)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, TPF = TPF)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, TPF = TPF)
            }

        resBoot <- simplify2array(resBoot)    
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        FPF.s <- simplify2array(resBoot["FPF.s",])
        if(length(TPF) == 1) {            	
        	thresholds.s <- matrix(thresholds.s, nrow = 1)
        	FPF.s <- matrix(FPF.s, nrow = 1)
        }

    } else {
        stop("nsave should be larger than zero.")
    }
    alpha <- (1-ci.level)/2
    np <- length(TPF)
    thresholds <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
	rownames(thresholds) <- TPF

	thresholds[,1] <- apply(thresholds.s, 1, mean)
	thresholds[,2] <- apply(thresholds.s, 1, quantile, prob = alpha)
	thresholds[,3] <- apply(thresholds.s, 1, quantile, prob = 1-alpha)

	FPF <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
	rownames(FPF) <- TPF

	FPF[,1] <- apply(FPF.s, 1, mean)
	FPF[,2] <- apply(FPF.s, 1, quantile, prob = alpha)
	FPF[,3] <- apply(FPF.s, 1, quantile, prob = 1-alpha)

	res <- list()
	res$thresholds <- thresholds
	res$FPF <- FPF
	res$TPF <- TPF
	res
}
