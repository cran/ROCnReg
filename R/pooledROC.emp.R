pooledROC.emp <-
function(marker, group, tag.h, data, p = seq(0,1,l=101), B = 1000, method = c("ncoutcome","coutcome"), pauc = pauccontrol(), parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    doBoostROC <- function(i, data, method, pauc, p) {
        data.boot <- bootstrap.sample(data, "group", method = method)
        yhb <- data.boot$y[data.boot$group == 0]
        ydb <- data.boot$y[data.boot$group == 1]
        obj.boot <- compute.ROC(yh = yhb, yd = ydb, pauc = pauc, p = p)

        res <- list()
        res$ROC <- obj.boot$ROC
        res$AUC  <- obj.boot$AUC
        if(pauc$compute) {
            res$pAUC  <- obj.boot$pAUC
        }
        res
    }

    compute.ROC <- function(yh, yd, pauc, p = seq(0,1,l=101)) {
    	n0 <- length(yh)
    	n1 <- length(yd)

        F1emp <- ecdf(yd)
        rocemp <- 1 - F1emp(quantile(yh, 1-p, type = 1))
        aucemp <- sum(outer(yh, yd, "<"))/(n0*n1) + sum(outer(yh, yd, "=="))/(2*n0*n1)
        
		if(pauc$compute) {
			if(pauc$focus == "FPF"){
				aux <- numeric(n1)
				for(j in 1:n1){
					aux[j] <- min(survival(yh, yd[j]), pauc$value)
				}  
				paucemp <- pauc$value - (1/n1)*sum(aux)
			} else {
				aux <- numeric(n0)
				for(i in 1:n0){
					aux[i] <- max(survival(yd, yh[i]), pauc$value)
				}
				paucemp <-(1/n0)*sum(aux) - pauc$value
			}
		}

        res <- list()
        res$p <- p
        res$ROC <- rocemp
        res$AUC <- aucemp
        if(pauc$compute) {
            res$pAUC <- paucemp
        }
        res
        
    }
    pauc <- do.call("pauccontrol", pauc)
    
    method <- match.arg(method)
    parallel <- match.arg(parallel)

    np <- length(p)

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
    
    res <- compute.ROC(yh = yh.wom, yd = yd.wom, pauc = pauc, p = p)
    rocemp <- res$ROC
    aucemp <- res$AUC
    paucemp <- res$pAUC
    
    data.original <- data.frame(y = c(yh.wom, yd.wom), group = c(rep(0, n0), rep(1, n1)))

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
                    parallel::mclapply(seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p)
            }

        resBoot <- simplify2array(resBoot)    
        rocempb <- simplify2array(resBoot["ROC",])
        aucempb <- unlist(resBoot["AUC",])
        if(pauc$compute){
            paucempb <- unlist(resBoot["pAUC",])
        }
    }
    columns <- switch(as.character(B>0),"TRUE" = 1:3,"FALSE"=1)
    col.names <- c("est","ql", "qh")[columns]
    
    poolROC <- matrix(0, ncol = length(columns), nrow = np, dimnames = list(1:np, col.names))
    poolROC[,1] <- rocemp
    AUC <- aucemp
    pAUC <- paucemp
    if(B > 0) {
        poolROC[,2] <- apply(rocempb, 1, quantile, prob = 0.025)
        poolROC[,3] <- apply(rocempb, 1, quantile, prob = 0.975)
        AUC <- c(AUC, quantile(aucempb,c(0.025,0.975)))
        if(pauc$compute) {
            pAUC <- c(pAUC, quantile(paucempb,c(0.025,0.975)))
        }
    }
    names(AUC) <- col.names
     if(pauc$compute){
        names(pAUC) <- col.names
    }
    
    res <- list()
    res$call <- match.call()
    res$marker <- list(h = yh, d = yd)
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$p <- p
    res$ROC <- poolROC
    res$AUC <- AUC
    if(pauc$compute) {
        res$pAUC <- if(pauc$focus == "FPF") {
        	pAUC/pauc$value
        } else {
        	pAUC/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    class(res) <- c("pooledROC.emp", "pooledROC")
    res
}
