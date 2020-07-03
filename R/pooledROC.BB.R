pooledROC.BB <-
function(marker, group, tag.h, data, p = seq(0, 1, l = 101), B = 5000, pauc = pauccontrol(), parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    compute.ROC <- function(i, yh, yd, pauc, p = seq(0,1,l=101)) {
        n0 <- length(yh)
        n1 <- length(yd)

        q <- rexp(n0, 1)
        weights.h <- q/sum(q)
        
        q1 <- rexp(n1, 1)
        weights.d <- q1/sum(q1)
        
        u <- 1 - ewcdf(yh, weights.h)(yd)
        rocbbpool <- ewcdf(u, weights.d)(p)
        aucbbpool <- sum((1-u)* weights.d)

        if(pauc$compute) {
            if(pauc$focus == "FPF") {
                paucbbpool <- sum((pauc$value - pmin(u, pauc$value))*weights.d)
            } else {
                u1 <- 1 - ewcdf(yd, weights.d)(yh)
                paucbbpool <- sum((pmax(u1, pauc$value)-pauc$value)*weights.h)
            }
        }
        res <- list()
        res$weights.d <- weights.d
        res$weights.h <- weights.h
        res$ROC <- rocbbpool
        res$AUC <- aucbbpool
        if(pauc$compute) {
            res$pAUC <- paucbbpool
        }
        res
    }

    pauc <- do.call("pauccontrol", pauc)
    
    parallel <- match.arg(parallel)
    
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
                    parallel::mclapply(seq_len(B), compute.ROC, yh = yh.wom, yd = yd.wom, pauc = pauc, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), compute.ROC, yh = yh.wom, yd = yd.wom, pauc = pauc, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), compute.ROC, yh = yh.wom, yd = yd.wom, pauc = pauc, p = p)
                        }                         
                    }
                }
            } else {
                lapply(seq_len(B), compute.ROC, yh = yh.wom, yd = yd.wom, pauc = pauc, p = p)
            }

        resBoot <- simplify2array(resBoot)
        weights.d <- simplify2array(resBoot["weights.d",])    
        weights.h <- simplify2array(resBoot["weights.h",])    
        rocbbpool <- simplify2array(resBoot["ROC",])
        aucbbpool <- unlist(resBoot["AUC",])
        if(pauc$compute){
            paucbbpool <- unlist(resBoot["pAUC",])
        }
    } else {
        stop("B should be larger than zero.")
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
