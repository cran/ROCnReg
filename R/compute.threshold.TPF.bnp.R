compute.threshold.TPF.bnp <-
function(object_h = NULL, object_d, newdata, TPF = 0.5, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    
    doMCMCTH <- function(k, object_h, object_d = NULL, Xhp = NULL, Xdp, TPF) {
        ncov <- nrow(Xhp)
        np <- length(TPF)
        
        # Thresholds
        thresholds <- matrix(0, ncol = ncov, nrow = np)
        if(is.null(object_d$probs)){
            for(incov in 1:ncov) {
                thresholds[,incov] <- qnorm(1-TPF, mean = Xdp[incov,]%*%object_d$beta[k,], sd = object_d$sd[k])            
            }
        } else {
            mu.d <- Xdp%*%t(object_d$beta[k,,])
            for(incov in 1:ncov) {
                aux <- norMix(mu = c(mu.d[incov,]), sigma = object_d$sd[k,], w = object_d$probs[k,])
                thresholds[,incov] <- qnorMix(1-TPF, aux)
            }
        }

        # FPF        
        if(!is.null(object_h)) {
            FPF <- matrix(0, ncol = ncov, nrow = np)
            ncov <- nrow(Xhp)
            if(is.null(object_h$probs)){
                mu.h <- Xhp%*%object_h$beta[k,]
                for(incov in 1:ncov) {
                     FPF[,incov] <- 1-pnorm(thresholds[,incov], mu.h[incov], object_h$sd[k])
                }
            } else {
                mu.h <- Xhp%*%t(object_h$beta[k,,])
                for(incov in 1:ncov) {
                    aux <- norMix(mu = c(mu.h[incov,]), sigma = object_h$sd[k,], w = object_h$probs[k,])
                    FPF[,incov] <- 1-pnorMix(thresholds[,incov], aux)
                }
            }
        }
        res <- list()
        res$thresholds <- thresholds
        if(!is.null(object_d)) {
            res$FPF <- FPF
        } else {
            res$foo <- 1
        }
        res
    }

    parallel <- match.arg(parallel)

    Xdp <- predict.design.matrix.bnp(object_d$mm, newdata)$X
    if(!is.null(object_h)) {
        Xhp <- predict.design.matrix.bnp(object_h$mm, newdata)$X
    } else {
        Xhp <- NULL
    }

    nrep <- nrow(object_d$beta)

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
                    parallel::mclapply(seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, TPF = TPF, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, TPF = TPF)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, TPF = TPF)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(nrep), doMCMCTH, object_h = object_h, object_d = object_d, Xhp = Xhp, Xdp = Xdp, TPF = TPF)
            }

        resBoot <- simplify2array(resBoot)
        thresholds <- simplify2array(resBoot["thresholds",])
        if(!is.null(object_d)) {
            FPF <- simplify2array(resBoot["FPF",])
        }

    } else {
        stop("nsave should be larger than zero.")
    }
    alpha <- (1-ci.level)/2
    np <- dim(thresholds)[1]
    ncov <- dim(thresholds)[2]
    
    thresholdsm <- thresholdsl <- thresholdsh <- matrix(0, nrow = np, ncol = ncov)
    rownames(thresholdsm) <- rownames(thresholdsl) <- rownames(thresholdsh) <- TPF
    
    thresholdsm <- apply(thresholds, c(1,2), mean)
    thresholdsl <- apply(thresholds, c(1,2), quantile, alpha)
    thresholdsh <- apply(thresholds, c(1,2), quantile, 1-alpha)
    
    # Organise results as desired
    thresholds.ret <- vector("list", length(TPF))
    names(thresholds.ret) <- TPF
    for(i in 1:length(TPF)){
        thresholds.ret[[i]] <- matrix(ncol = 3, nrow = ncov)
        thresholds.ret[[i]][,1] <- thresholdsm[i,] 
        thresholds.ret[[i]][,2] <- thresholdsl[i,]
        thresholds.ret[[i]][,3] <- thresholdsh[i,]
        colnames(thresholds.ret[[i]]) <- c("est", "ql", "qh")
    }
    res <- list()
    res$thresholds <- thresholds.ret

    if(!is.null(object_h)) {
        # Organise results as desired
        FPFm <- apply(FPF, c(1,2), mean)
        FPFl <- apply(FPF, c(1,2), quantile, alpha)
        FPFh <- apply(FPF, c(1,2), quantile, 1-alpha)

        # Organised results as desired
        FPF.ret <- vector("list", length(TPF))
        names(FPF.ret) <- TPF
        for(i in 1:length(TPF)){
            FPF.ret[[i]] <- matrix(ncol = 3, nrow = ncov)

            FPF.ret[[i]][,1] <- FPFm[i,] 
            FPF.ret[[i]][,2] <- FPFl[i,]
            FPF.ret[[i]][,3] <- FPFh[i,]

            colnames(FPF.ret[[i]]) <- c("est", "ql", "qh")
        }
        res$FPF <- FPF.ret
    }
    res
}
