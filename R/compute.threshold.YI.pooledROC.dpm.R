compute.threshold.YI.pooledROC.dpm <-
function(object, ci.level = 0.95, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    
    doMCMCTH <- function(k, res0, res1, grid) {
        p0 <- res0$probs
        p1 <- res1$probs

        if(is.null(p0) & is.null(p1)){
            F0 <- pnorm(grid, mean = res0$mu[k], sd = res0$sd[k])
            F1 <- pnorm(grid, mean = res1$mu[k], sd = res1$sd[k])
            
            dif <- F0 - F1
            
            thresholds.s <- mean(grid[which(dif == max(dif))])
            YI.s <- max(dif)
            
            TPF.s <- 1 - pnorm(thresholds.s, mean = res1$mu[k], sd = res1$sd[k])
            FPF.s <- 1 - pnorm(thresholds.s, mean = res0$mu[k], sd = res0$sd[k])
        } else if(is.null(p0) & !is.null(p1)){
            aux1 <- norMix(mu = res1$mu[k,], sigma = res1$sd[k,], w = p1[k,])
            
            F0 <- pnorm(grid, mean = res0$mu[k], sd = res0$sd[k])
            F1 <- pnorMix(grid, aux1)
            
            dif <- F0 - F1
            
            thresholds.s <- mean(grid[which(dif == max(dif))])
            YI.s <- max(dif)
            
            TPF.s <- 1 - pnorMix(thresholds.s, aux1)
            FPF.s <- 1 - pnorm(thresholds.s, mean = res0$mu[k], sd = res0$sd[k])
            
        } else if (!is.null(p0) & is.null(p1)){
            aux0 <- norMix(mu = res0$mu[k,], sigma = res0$sd[k,], w = p0[k,])
            
            F1 <- pnorm(grid, mean = res1$mu[k], sd = res1$sd[k])
            F0 <- pnorMix(grid, aux0)
            
            dif <- F0 - F1
            
            thresholds.s <- mean(grid[which(dif == max(dif))])
            YI.s <- max(dif)
            
            TPF.s <- 1 - pnorm(thresholds.s, mean = res1$mu[k], sd = res1$sd[k])
            FPF.s <- 1 - pnorMix(thresholds.s, aux0)
        } else {
            aux0 <- norMix(mu = res0$mu[k,], sigma = res0$sd[k,], w = p0[k,])
            aux1 <- norMix(mu = res1$mu[k,], sigma = res1$sd[k,], w = p1[k,])
            
            F0 <- pnorMix(grid, aux0)
            F1 <- pnorMix(grid, aux1)
            
            dif <- F0 - F1
            
            thresholds.s <- mean(grid[which(dif == max(dif))])
            YI.s <- max(dif)
            TPF.s <- 1 - pnorMix(thresholds.s, aux1)
            FPF.s <- 1 - pnorMix(thresholds.s, aux0)
        }

        res <- list()
        res$thresholds.s <- thresholds.s
        res$YI.s <- YI.s
        res$TPF.s <- TPF.s
        res$FPF.s <- FPF.s
        res

    }

    if(class(object)[1] != "pooledROC.dpm") {
        stop(paste0("This function can not be used for this object class: ", class(object)[1]))
    }

    parallel <- match.arg(parallel)

    y0 <- object$marker$h[!object$missing.ind$h]
    y1 <- object$marker$d[!object$missing.ind$d]

    grid <- seq(min(c(y0,y1))-1, max(c(y0,y1))+1, len = max(c(length(y0), length(y1)), 500))
    ngrid <- length(grid)

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
                    parallel::mclapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, grid = grid, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, grid = grid)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, grid = grid)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, grid = grid)
            }

        resBoot <- simplify2array(resBoot)    
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        YI.s <- simplify2array(resBoot["YI.s",])
        TPF.s <- simplify2array(resBoot["TPF.s",])
        FPF.s <- simplify2array(resBoot["FPF.s",])

    } else {
        stop("nsave should be larger than zero.")
    }

    alpha <- (1-ci.level)/2

    thresholds <- c(mean(thresholds.s), quantile(thresholds.s, c(alpha, 1-alpha)))
    YI <- c(mean(YI.s), quantile(YI.s, c(alpha, 1-alpha)))
    FPF <- c(mean(FPF.s), quantile(FPF.s, c(alpha, 1-alpha)))
    TPF <- c(mean(TPF.s), quantile(TPF.s, c(alpha, 1-alpha)))
    
    names(thresholds) <- names(YI) <- names(FPF) <- names(TPF) <- c("est","ql", "qh")
    
    res <- list()
    res$call <- match.call()
    res$thresholds <- thresholds
    res$YI <- YI
    res$FPF <- FPF
    res$TPF <- TPF
    res
}
