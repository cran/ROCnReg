compute.threshold.FPF.bnp <-
function(object_h, object_d = NULL, newdata, FPF = 0.5, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    
    doMCMCTH <- function(k, object_h, object_d = NULL, Xhp, Xdp = NULL, FPF) {
        ncov <- nrow(Xhp)
        np <- length(FPF)
        
        # Thresholds
        thresholds <- matrix(0, ncol = ncov, nrow = np)
        if(is.null(object_h$probs)){
            for(incov in 1:ncov) {
                thresholds[,incov] <- qnorm(1-FPF, mean = Xhp[incov,]%*%object_h$beta[k,], sd = object_h$sd[k])            
            }
        } else {
            mu.h <- Xhp%*%t(object_h$beta[k,,])
            for(incov in 1:ncov) {
                aux <- norMix(mu = c(mu.h[incov,]), sigma = object_h$sd[k,], w = object_h$probs[k,])
                thresholds[,incov] <- qnorMix(1-FPF, aux)
            }
        }

        # TPF        
        if(!is.null(object_d)) {
            TPF <- matrix(0, ncol = ncov, nrow = np)
            ncov <- nrow(Xdp)            
            if(is.null(object_d$probs)){
                mu.d <- Xdp%*%object_d$beta[k,]
                for(incov in 1:ncov) {
                     TPF[,incov] <- 1-pnorm(thresholds[,incov], mu.d[incov], object_d$sd[k])            
                }
            } else {
                mu.d <- Xdp%*%t(object_d$beta[k,,])
                for(incov in 1:ncov) {
                    aux <- norMix(mu = c(mu.d[incov,]), sigma = object_d$sd[k,], w = object_d$probs[k,])
                    TPF[,incov] <- 1-pnorMix(thresholds[,incov], aux)
                }
            }
        }
        res <- list()
        res$thresholds <- thresholds
        if(!is.null(object_d)) {
            res$TPF <- TPF
        } else {
            res$foo <- 1
        }
        res
    }

    parallel <- match.arg(parallel)

    Xhp <- predict.design.matrix.bnp(object_h$mm, newdata)$X
    if(!is.null(object_d)) {
        Xdp <- predict.design.matrix.bnp(object_d$mm, newdata)$X
    } else {
        Xdp <- NULL
    }

    nrep <- nrow(object_h$beta)

    if(nrep > 0) {
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
                    parallel::mclapply(seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, FPF = FPF, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, FPF = FPF)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, FPF = FPF)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, FPF = FPF)
            }

        resBoot <- simplify2array(resBoot)
        thresholds <- simplify2array(resBoot["thresholds",])
        if(!is.null(object_d)) {
            TPF <- simplify2array(resBoot["TPF",])
        }

    } else {
        stop("nsave should be larger than zero.")
    }

    np <- dim(thresholds)[1]
    ncov <- dim(thresholds)[2]
    
    thresholdsm <- thresholdsl <- thresholdsh <- matrix(0, nrow = np, ncol = ncov)
    rownames(thresholdsm) <- rownames(thresholdsl) <- rownames(thresholdsh) <- FPF
    
    thresholdsm <- apply(thresholds, c(1,2), mean)
    thresholdsl <- apply(thresholds, c(1,2), quantile, 0.025)
    thresholdsh <- apply(thresholds, c(1,2), quantile, 0.975)
    
    # Organise results as desired
    thresholds.ret <- vector("list", length(FPF))
    names(thresholds.ret) <- FPF
    for(i in 1:length(FPF)){
        thresholds.ret[[i]] <- matrix(ncol = 3, nrow = ncov)
        thresholds.ret[[i]][,1] <- thresholdsm[i,] 
        thresholds.ret[[i]][,2] <- thresholdsl[i,]
        thresholds.ret[[i]][,3] <- thresholdsh[i,]
        colnames(thresholds.ret[[i]]) <- c("est", "ql", "qh")
    }
    res <- list()
    res$thresholds <- thresholds.ret

    if(!is.null(object_d)) {
        # Organise results as desired
        TPFm <- apply(TPF, c(1,2), mean)
        TPFl <- apply(TPF, c(1,2), quantile, 0.025)
        TPFh <- apply(TPF, c(1,2), quantile, 0.975)

        # Organised results as desired
        TPF.ret <- vector("list", length(FPF))
        names(TPF.ret) <- FPF
        for(i in 1:length(FPF)){
            TPF.ret[[i]] <- matrix(ncol = 3, nrow = ncov)

            TPF.ret[[i]][,1] <- TPFm[i,] 
            TPF.ret[[i]][,2] <- TPFl[i,]
            TPF.ret[[i]][,3] <- TPFh[i,]

            colnames(TPF.ret[[i]]) <- c("est", "ql", "qh")
        }
        res$TPF <- TPF.ret
    }
    res
}
