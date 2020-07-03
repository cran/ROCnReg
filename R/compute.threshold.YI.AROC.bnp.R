compute.threshold.YI.AROC.bnp <-
function(object, newdata, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    doMCMCTH <- function(k, res0, L, yd, Xd, Xhp, p) {
        nd <- length(yd)
        np <- length(p)
        npred <- nrow(Xhp)

        if(L == 1) {
            up <- 1 - pnorm(yd, mean = Xd%*%res0$beta[k,], sd = res0$sd[k])
        } else {
            up <- 1 - apply(t(res0$probs[k,]*t(pnorm(yd, mean = tcrossprod(Xd, res0$beta[k,,]), sd = rep(res0$sd[k,], each = length(yd))))), 1, sum)
        }

        aux <- rexp(nd, 1)
        weights <- aux/sum(aux)

        AROC <- apply(weights*outer(up, p, "<="), 2, sum)

        dif <- AROC -  p
        FPF.s <- mean(p[which(dif == max(dif))])
        YI.s <- max(dif)
        
        thresholds <- vector(length = npred)

        if(L > 1){
            mu.h <- Xhp%*%t(res0$beta[k,,])
            for(i in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[i,]), sigma = res0$sd[k,], w = res0$probs[k,])
                thresholds[i] <- qnorMix(1 - FPF.s, aux0)
            }
        }
        if(L == 1){
            mu.h <- Xhp%*%res0$beta[k,]
            #for(i in 1:npred) {
            #    thresholds[i] <- qnorm(1 - FPF.s, mu.h[i], Sigma0[k])
            #}
            thresholds <- qnorm(1 - FPF.s, mu.h, res0$sd[k])
        }

        res <- list()
        res$FPF.s <- FPF.s
        res$YI.s <- YI.s
        res$thresholds <- thresholds
        res
    }


    if(class(object)[1] != "AROC.bnp") {
        stop(paste0("This function cannot be used for this object class: ", class(object)[1]))
    }

    parallel <- match.arg(parallel)
    
    #names.cov <- all.vars(object$fit$formula)[-1]
    names.cov <- get_vars_formula(object$fit$formula)
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    
    if(missing(newdata)) {
        newdata <- cROCData(object$data, names.cov, object$group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop = FALSE])
    }

    p <- seq(0, 1, length = 500)

    X0p <- predict(object$fit$mm, newdata = newdata)$X

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
                    parallel::mclapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit, L = object$prior$L, yd = object$data_model$y$d, Xd = object$data_model$X$d, Xhp = X0p, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit, L = object$prior$L, yd = object$data_model$y$d, Xd = object$data_model$X$d, Xhp = X0p, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit, L = object$prior$L, yd = object$data_model$y$d, Xd = object$data_model$X$d, Xhp = X0p, p = p)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(object$mcmc$nsave), doMCMCTH, res0 = object$fit, L = object$prior$L, yd = object$data_model$y$d, Xd = object$data_model$X$d, Xhp = X0p, p = p)
            }

        resBoot <- simplify2array(resBoot)    
        FPF.s <- simplify2array(resBoot["FPF.s",])
        YI.s <- simplify2array(resBoot["YI.s",])
        thresholds <- simplify2array(resBoot["thresholds",])

    } else {
        stop("nsave should be larger than zero.")
    }
    
    YI <- c(mean(YI.s), quantile(YI.s, c(0.025,0.975)))
    FPF <- c(mean(FPF.s), quantile(FPF.s, c(0.025,0.975)))
    names(YI) <- names(FPF) <- c("est","ql", "qh")
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$thresholds <- cbind(est = apply(thresholds, 1, mean), ql = apply(thresholds, 1, quantile, 0.025), qh = apply(thresholds, 1, quantile, 0.975))
    res$YI <- YI
    res$FPF <- FPF
    res
}
