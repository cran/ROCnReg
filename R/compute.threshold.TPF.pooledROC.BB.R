compute.threshold.TPF.pooledROC.BB <-
function(object, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
	doMCMCTH <- function(k, y.h, weights.h, y.d, weights.d, TPF) {
		thresholds.s <- quantile(spatstat.univar::ewcdf(y.d, weights.d[,k]), 1- TPF, type = 1)
		FPF.s <- 1 - spatstat.univar::ewcdf(y.h, weights.h[,k])(thresholds.s)
		res <- list()
		res$thresholds.s <- thresholds.s
		res$FPF.s <- FPF.s
		res
	}

	if(class(object)[1] != "pooledROC.BB") {
		stop(paste0("This function can not be used for this object class: ", class(object)[1]))
	}

	parallel <- match.arg(parallel)

	B <- ncol(object$weights$h)

	y.h <- object$marker$h[!object$missing.ind$h]
	y.d <- object$marker$d[!object$missing.ind$d]

	weights.h <- object$weights$h
	weights.d <- object$weights$d

	if(B > 0) {        
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
                    parallel::mclapply(seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, TPF = TPF , mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, TPF = TPF )
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, TPF = TPF)
                        }                         
                    }
                }
            } else {
                lapply(seq_len(B), doMCMCTH, y.h = y.h, weights.h = weights.h, y.d = y.d, weights.d = weights.d, TPF = TPF)
            }

        resBoot <- simplify2array(resBoot)
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        FPF.s <- simplify2array(resBoot["FPF.s",])
        if(length(FPF) == 1) {            	
        	thresholds.s <- matrix(thresholds.s, nrow = 1)
        	FPF.s <- matrix(FPF.s, nrow = 1)
        }
    } else {
        stop("B should be larger than zero.")
    }
    alpha <- (1-ci.level)/2
	np <- length(FPF)
	thresholds <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
	rownames(thresholds) <- TPF

	thresholds[,1] <- apply(thresholds.s, 1, mean)
	thresholds[,2] <- apply(thresholds.s, 1, quantile, prob = alpha)
	thresholds[,3] <- apply(thresholds.s, 1, quantile, prob = 1 - alpha)

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
