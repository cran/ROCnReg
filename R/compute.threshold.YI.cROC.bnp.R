compute.threshold.YI.cROC.bnp <-
function(object, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    
    doMCMCTH <- function(k, res0, res1, Xhp, Xdp, Lh, Ld, grid) {
        npred <- nrow(Xhp)
        ngrid <- length(grid)

        F0 <- F1 <- matrix(0, nrow = ngrid, ncol = npred)
        thresholds.s <- YI.s <- TPF.s <- FPF.s <- vector(length = npred)

        if(Ld == 1 & Lh == 1){
            mu.h <- Xhp%*%res0$beta[k,]
            mu.d <- Xdp%*%res1$beta[k,]
            
            for(l in 1:npred){
                F0[,l] <- pnorm(grid, mu.h[l], res0$sd[k])
                F1[,l] <- pnorm(grid, mu.d[l], res1$sd[k])
                
                dif <- abs(F0[,l] -  F1[,l])
                thresholds.s[l] <- mean(grid[which(dif == max(dif))])
                YI.s[l] <- max(dif)
                
                TPF.s[l] <- 1 - pnorm(thresholds.s[l], mu.d[l], res1$sd[k])
                FPF.s[l] <- 1 - pnorm(thresholds.s[l], mu.h[l], res0$sd[k])
            }
        }
        
        if(Ld == 1 & Lh > 1){
            mu.d <- Xdp%*%res1$beta[k,]
            mu.h <- Xhp%*%t(res0$beta[k,,])
            for(l in 1:npred){
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = res0$sd[k,], w = res0$probs[k,])
                F0[,l] <- pnorMix(grid, aux0)
                F1[,l] <- pnorm(grid, mu.d[l], res1$sd[k])
                
                dif <- abs(F0[,l] -  F1[,l])
                thresholds.s[l] <- mean(grid[which(dif == max(dif))])
                YI.s[l] <- max(dif)
                
                TPF.s[l] <- 1 - pnorm(thresholds.s[l], mu.d[l], res1$sd[k])
                FPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux0)
            }
        }
        
        if(Ld > 1 & Lh == 1){
            mu.h <- Xhp%*%res0$beta[k,]
            mu.d <- Xdp%*%t(res1$beta[k,,])
            for(l in 1:npred){
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = res1$sd[k,], w = res1$probs[k,])
                F0[,l] <- pnorm(grid, mu.h[l], res0$sd[k])
                F1[,l] <- pnorMix(grid, aux1)
                
                dif <- abs(F0[,l] -  F1[,l])
                thresholds.s[l] <- mean(grid[which(dif == max(dif))])
                YI.s[l] <- max(dif)
                
                TPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux1)
                FPF.s[l] <- 1 - pnorm(thresholds.s[l], mu.h[l], res0$sd[k])
            }
        }
        
        if(Ld > 1 & Lh > 1){
            mu.h <- Xhp%*%t(res0$beta[k,,])
            mu.d <- Xdp%*%t(res1$beta[k,,])
            for(l in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = res0$sd[k,], w = res0$probs[k,])
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = res1$sd[k,], w = res1$probs[k,])
                
                F0[,l] <- pnorMix(grid, aux0)
                F1[,l] <- pnorMix(grid, aux1)
                
                difbb <- abs(F0[,l] -  F1[,l])
                thresholds.s[l] <- mean(grid[which(difbb == max(difbb))])
                YI.s[l] <- max(difbb)
                
                TPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux1)
                FPF.s[l] <- 1 - pnorMix(thresholds.s[l], aux0)
            }
        }

        res <- list()
        res$thresholds.s <- thresholds.s
        res$YI.s <- YI.s
        res$TPF.s <- TPF.s
        res$FPF.s <- FPF.s
        res

    }

    if(class(object)[1] != "cROC.bnp") {
        stop(paste0("This function can not be used for this object class: ", class(object)[1]))
    }

    parallel <- match.arg(parallel)

    #names.cov.h <- all.vars(object$fit$h$formula)[-1]
    #names.cov.d <- all.vars(object$fit$d$formula)[-1]
    #names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])

    names.cov.h <- get_vars_formula(object$fit$h$formula)
    names.cov.d <- get_vars_formula(object$fit$d$formula)
    names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])
    
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    
    if(missing(newdata)) {
        newdata <- cROCData(object$data, names.cov, object$group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop=FALSE])
    }
    
    # Compute F_D|X and F_{\bar{D}}|X
    X0p <- predict(object$fit$h$mm, newdata = newdata)$X
    X1p <- predict(object$fit$d$mm, newdata = newdata)$X
    
    Lh <- object$prior$h$L
    Ld <- object$prior$d$L
    
    y0 <- object$data_model$y$h
    y1 <- object$data_model$y$d
    n0 <- length(y0)
    n1 <- length(y1)
    
    grid  <- seq(min(c(y0, y1), na.rm = TRUE) - 1, max(c(y0, y1), na.rm = TRUE) + 1, length = max(500, c(n0,n1)))


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
                    parallel::mclapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, Lh = Lh, Ld = Ld, Xhp = X0p, Xdp = X1p, grid = grid, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, Lh = Lh, Ld = Ld, Xhp = X0p, Xdp = X1p, grid = grid)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, Lh = Lh, Ld = Ld, Xhp = X0p, Xdp = X1p, grid = grid)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit$h, res1 = object$fit$d, Lh = Lh, Ld = Ld, Xhp = X0p, Xdp = X1p, grid = grid)
            }
            
        resBoot <- simplify2array(resBoot)    
        thresholds.s <- simplify2array(resBoot["thresholds.s",])
        YI.s <- simplify2array(resBoot["YI.s",])
        TPF.s <- simplify2array(resBoot["TPF.s",])
        FPF.s <- simplify2array(resBoot["FPF.s",])

    } else {
        stop("nsave should be larger than zero.")
    }

    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$thresholds <- cbind(est = apply(thresholds.s, 1, mean), ql = apply(thresholds.s, 1, quantile, 0.025), qh = apply(thresholds.s, 1, quantile, 0.975))
    res$YI <- cbind(est = apply(YI.s, 1, mean), ql = apply(YI.s, 1, quantile, 0.025), qh = apply(YI.s, 1, quantile, 0.975))
    res$FPF <- cbind(est = apply(FPF.s, 1, mean), ql = apply(FPF.s, 1, quantile, 0.025), qh = apply(FPF.s, 1, quantile, 0.975))
    res$TPF <- cbind(est = apply(TPF.s, 1, mean), ql = apply(TPF.s, 1, quantile, 0.025), qh = apply(TPF.s, 1, quantile, 0.975))
    res
}
