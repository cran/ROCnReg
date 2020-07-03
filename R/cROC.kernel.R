cROC.kernel <-
function(marker, covariate, group, tag.h, bw = c("LS","AIC"), regtype = c("LC","LL"), data, newdata, pauc = pauccontrol(), p = seq(0,1,l = 101), B = 1000, parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    doBoostROC <- function(i, marker, covariate, group, tag.h, bw, regtype, croc, newdata, pauc, p) {
        
        data.boot.d <- croc$data.d
        data.boot.h <- croc$data.h

        res.h.b <- sample(croc$fit$h$fit.mean$resid/sqrt(croc$fit$h$fit.var$mean), replace = TRUE)
        res.d.b <- sample(croc$fit$d$fit.mean$resid/sqrt(croc$fit$d$fit.var$mean), replace = TRUE)

        data.boot.h[,marker] <- croc$fit$h$fit.mean$mean + sqrt(croc$fit$h$fit.var$mean)*res.h.b
        data.boot.d[,marker] <- croc$fit$d$fit.mean$mean + sqrt(croc$fit$d$fit.var$mean)*res.d.b
        data.boot <- rbind(data.boot.d, data.boot.h)

        obj.boot <- compute.ROC(marker, covariate, group, tag.h, bw, regtype, data.boot, newdata, pauc, p)
        
        res <- list()
        res$ROC <- obj.boot$ROC
        res$AUC <- obj.boot$AUC
        if(pauc$compute){
            res$pAUC <- obj.boot$pAUC
        }
        res
    }

    compute.ROC <- function(marker, covariate, group, tag.h, bw, regtype, data, newdata, pauc, p = seq(0,1,l = 101)) {
        data.h <- data[data[,group] == tag.h,]
        data.d <- data[data[,group] != tag.h,]
        
        n0 <- nrow(data.h)
        n1 <- nrow(data.d)
        
        np <- length(p)
        
        x0 <- data.h[,covariate]
        y0 <- data.h[,marker]
        x1 <- data.d[,covariate]
        y1 <- data.d[,marker]
        
        xp <- newdata[, covariate]
        
        # Fit the location-scale model in the healthy population
        bw.mean.h <-  npregbw(ydat = y0, xdat = x0, regtype = regtype, bwmethod = bw)
        fit.mean.h <- npreg(bw.mean.h, exdat = x0,residuals = TRUE)
        bw.var.h <- npregbw(ydat = (fit.mean.h$resid^2), xdat = x0, regtype = "lc", bwmethod = bw)
        fit.var.h <- npreg(bw.var.h, exdat = x0, residuals = TRUE)
        
        res0p <- fit.mean.h$resid/sqrt(fit.var.h$mean)
        F0res <- ecdf(res0p)
        
        # Fit the location-scale model in the diseased population
        bw.mean.d <-  npregbw(ydat = y1, xdat = x1, regtype = regtype, bwmethod = bw)
        fit.mean.d <- npreg(bw.mean.d, exdat = x1, residuals = TRUE)
        bw.var.d <- npregbw(ydat = (fit.mean.d$resid^2), xdat = x1, regtype = "lc", bwmethod = bw)
        fit.var.d <- npreg(bw.var.d, exdat = x1, residuals = TRUE)
        
        res1p <- fit.mean.d$resid/sqrt(fit.var.d$mean)
        F1res <- ecdf(res1p)
        
        # Evaluate the model in the newdata, and compute the cROC
        fit.mean.h.p <- npreg(bw.mean.h, exdat = xp,residuals = TRUE)
        fit.var.h.p <- npreg(bw.var.h, exdat = xp, residuals = TRUE)
        
        fit.mean.d.p <- npreg(bw.mean.d, exdat = xp,residuals = TRUE)
        fit.var.d.p <- npreg(bw.var.d, exdat = xp, residuals = TRUE)
        # Three terms of the cROC
        a <- (fit.mean.h.p$mean - fit.mean.d.p$mean)/sqrt(fit.var.d.p$mean)
        b <- sqrt(fit.var.h.p$mean)/sqrt(fit.var.d.p$mean)
        
        #cdf0 <- apply(outer(res0p, res0p, "<="), 2, mean)
        #cdf0_inv <- apply(outer(cdf0, 1-p, ">="), 2, function(x, z) {
        #        res <- min(z[x]) #min(c(z[x], max(z)))
        #        res
        #    }, z = res0p)
        #cdf0_inv <- replace(cdf0_inv, is.infinite(cdf0_inv), max(res0p))
        cdf0_inv <- quantile(res0p, 1-p, type = 1)
        
        # cROC: rows covariates / columns FPF
        cROC <- 1 - apply(a + outer(b, cdf0_inv, "*"), c(1,2), F1res)
        cAUC <- apply(cROC, 1, simpson, set.p = p)

		if(pauc$compute){
			if(pauc$focus=="FPF"){
				pu <- seq(0, pauc$value, len = np)
				cdf0_inv_pauc <- quantile(res0p, 1-pu, type = 1)
				cROC_pauc <- 1 - apply(a + outer(b, cdf0_inv_pauc, "*"), c(1,2), F1res)
				cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
			} else{
				a <- (fit.mean.d.p$mean - fit.mean.h.p$mean)/sqrt(fit.var.h.p$mean)
				b <- sqrt(fit.var.d.p$mean)/sqrt(fit.var.h.p$mean)
				pu <- seq(pauc$value, 1, len = np)
				cdf1_inv_pauc <- quantile(res1p, 1-pu, type = 1)
				cROC_pauc <- apply(a + outer(b, cdf1_inv_pauc, "*"), c(1,2), F0res)
				cpAUC <- apply(cROC_pauc, 1, simpson, set.p = pu)
			}
		}

        res <- list()
        res$p <- p
        res$ROC <- cROC
        res$AUC <- cAUC
        if(pauc$compute){
            res$pAUC <- cpAUC
        }
        res$data.h <- data.h
        res$data.d <- data.d
        res$fit <- list()
        res$fit$h <- list(bw.mean = bw.mean.h, bw.var = bw.var.h, fit.mean = fit.mean.h, fit.var = fit.var.h)
        res$fit$d <- list(bw.mean = bw.mean.d, bw.var = bw.var.d, fit.mean = fit.mean.d, fit.var = fit.var.d)
        res
    }

    np <- length(p)
    # Check if the seq of FPF is correct
    if(min(p) != 0 | max(p) != 1 | np%%2 == 0) {
        warning("The set of FPFs at which the covariate-specific ROC curve is estimated is not correct. The set is used for calculating the AUC using Simpson's rule. As such, it should (a) start in 0; (b) end in 1; and (c) have an odd length.")
    }
    
    pauc <- do.call("pauccontrol", pauc)

    parallel <- match.arg(parallel)
    bw <- match.arg(bw)
    regtype <- match.arg(regtype)
    
    bw.aux <- switch(bw, "LS" = "cv.ls", "AIC" = "cv.aic")
    regtype.aux <- switch(regtype, "LC" = "lc", "LL" = "ll")
    
    names.cov <- covariate
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(sum(is.na(match(c(marker,names.cov,group), names(data)))))
        stop("Not all needed variables are supplied in data")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    if(length(unique(data[,group]))!= 2)
        stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep="")

    # New data, removing missing values
    data.new <- data[,c(marker,group,covariate)]
    omit.h <- apply(data.new[data.new[,group] == tag.h, c(marker, group, covariate)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.h, c(marker, group, covariate)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.h,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.h,,drop = FALSE][!omit.d,,drop = FALSE])
    
    # New data (for predictions)    
    if(missing(newdata)) {
        newdata <- cROCData(data.new, covariate, group)
    } else {
        newdata <- na.omit(newdata[,covariate,drop = FALSE])
    }
    
    croc <- compute.ROC(marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, data = data.new, newdata = newdata, pauc = pauc, p = p)
    crocp <- croc$ROC
    caucp <- croc$AUC
    if(pauc$compute){
        cpaucp <- croc$pAUC
    }
    
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
                    parallel::mclapply(seq_len(B), doBoostROC, marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, croc = croc, newdata = newdata, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), doBoostROC, marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, croc = croc, newdata = newdata, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), doBoostROC, marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, croc = croc, newdata = newdata, p = p)
                        }                          
                    }
                }
            } else {
                lapply(seq_len(B), doBoostROC, marker = marker, covariate = covariate, group = group, tag.h = tag.h, bw = bw.aux, regtype = regtype.aux, pauc = pauc, croc = croc, newdata = newdata, p = p)
            }

        resBoot <- simplify2array(resBoot)    
        crocpb <- simplify2array(resBoot["ROC",])
        caucb <- simplify2array(resBoot["AUC",])
        if(pauc$compute){
            cpaucb <- simplify2array(resBoot["pAUC",])
        }
    }
    columns <-switch(as.character(B > 0),"TRUE" = 1:3, "FALSE" = 1)
    col.names <-c("AUC","AUCql", "AUCqh")[columns]
    if(pauc$compute){
        col.names_pAUC <-c("pAUC","pAUCql", "pAUCqh")[columns]
    }
    cROCres <- list()
    cROCres$est <- crocp
    cAUCres <- matrix(0, ncol = length(columns), nrow = nrow(newdata), dimnames = list(1:nrow(newdata), col.names))
    cAUCres[,1] <- caucp
    if(pauc$compute){
        cpAUCres <- matrix(0, ncol = length(columns), nrow = nrow(newdata), dimnames = list(1:nrow(newdata), col.names))
        cpAUCres[,1] <- cpaucp
    }
    if(B > 0) {
        cROCres$ql <- apply(crocpb, c(1,2),quantile, prob=0.025)
        cROCres$qh <- apply(crocpb, c(1,2),quantile, prob=0.975)
        
        cAUCres[,2] <- apply(caucb, 1, quantile, prob=0.025)
        cAUCres[,3] <- apply(caucb, 1, quantile, prob=0.975)
        
        if(pauc$compute){
            cpAUCres[,2] <- apply(cpaucb, 1, quantile, prob=0.025)
            cpAUCres[,3] <- apply(cpaucb, 1, quantile, prob=0.975)
        }
    }
    colnames(cAUCres) <- col.names
    if(pauc$compute){
        colnames(cpAUCres) <- col.names_pAUC
    }
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$data <- data
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$marker <- marker
    res$group <- group
    res$tag.h <- tag.h
    res$covariate <- covariate
    res$p <- p
    res$ROC <- cROCres
    res$AUC <- cAUCres
    if(pauc$compute){
        res$pAUC <- if(pauc$focus == "FPF") {
        	cpAUCres/pauc$value
        } else {
        	cpAUCres/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    res$fit <- croc$fit
    res$ci.fit <- ifelse(B, TRUE, FALSE)
    class(res) <- c("cROC.kernel", "cROC")
    res
}
