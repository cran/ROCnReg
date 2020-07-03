pooledROC.kernel <-
function(marker, group, tag.h, data, p = seq(0, 1, l = 101),  bw = c("SRT","UCV"), B = 1000, method = c("ncoutcome","coutcome"), pauc = pauccontrol(), parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL) {
    
    doBoostROC <- function(i, data, method, pauc, p, bw) {
        data.boot <- bootstrap.sample(data, "group", method = method)
        yhb <- data.boot$y[data.boot$group == 0]
        ydb <- data.boot$y[data.boot$group == 1]
        obj.boot <- compute.ROC(yh = yhb, yd = ydb, pauc = pauc, p = p, bw = bw)

        res <- list()
        res$ROC <- obj.boot$ROC
        res$AUC  <- obj.boot$AUC
        if(pauc$compute) {
            res$pAUC  <- obj.boot$pAUC
        }
        res
    }
    compute.ROC <- function(yh, yd, pauc, p, bw) {
        if(bw == "SRT") {
            h0 <- suppressWarnings(bw.nrd0(x = yh))
            h1 <- suppressWarnings(bw.nrd0(x = yd))
        } else if (bw == "UCV") {
            h0 <- suppressWarnings(bw.ucv(x = yh))
            h1 <- suppressWarnings(bw.ucv(x = yd))
        } else {
            stop("Unknown bandwidth selector")
        }
        
        np <- length(p)
        rock <- numeric(np)
        for(j in 1:np){
            rock[j] <- 1-Gk(qFk(1-p[j], y = yh, h = h0), y = yd, h = h1)
        }
        auck <- simpson(ROC = rock, set.p = p)
		if(pauc$compute){
			if(pauc$focus == "FPF") {
				rock_pauc <- numeric(np)
		      	pu <- seq(0, pauc$value, len = np)
		      	for(j in 1:np){
		        	rock_pauc[j] <- 1-Gk(qFk(1-pu[j], y = yh, h = h0), y = yd, h = h1) 
		      	}  
		      	pauck <- simpson(ROC = rock_pauc, set.p = pu)
	      	} else {
	        	rock_pauc <- numeric(np)
	        	pu <- seq(pauc$value, 1, len = np)
	        	for(j in 1:np){
	          		rock_pauc[j] <- Gk(qFk(1-pu[j], y = yd, h = h1), y = yh, h = h0) 
	        	} 
	        	pauck <- simpson(ROC = rock_pauc, set.p = pu)
	    	}
	    }        
        res <- list()
        res$bws <- list(h = h0, d = h1)
        res$p <- p
        res$ROC <- rock
        res$AUC <- auck
        if(pauc$compute){
            res$pAUC <- pauck
        }
        res
    }
    pauc <- do.call("pauccontrol", pauc)

    parallel <- match.arg(parallel)
    method <- match.arg(method)
    bw <- match.arg(bw)

    np <- length(p)

    # Check if the seq of FPF is correct
    if(min(p) != 0 | max(p) != 1 | np%%2 == 0) {
        warning("The set of FPFs at which the pooled ROC curve is estimated is not correct. The set is used for calculating the AUC using Simpson's rule. As such, it should (a) start in 0; (b) end in 1; and (c) have an odd length.")
    }
    
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
    
    obj <- compute.ROC(yh = yh.wom, yd = yd.wom, pauc = pauc, p = p, bw = bw)
    rock <- obj$ROC
    auck <- obj$AUC
    if(pauc$compute){
        pauck <- obj$pAUC
    }
    data.original <- data.frame(y = c(yh.wom, yd.wom), group = c(rep(0,n0), rep(1,n1)))

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
                    parallel::mclapply(seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p, bw = bw, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p, bw = bw)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p, bw = bw)
                        }
                    }
                }
            } else {
                lapply(seq_len(B), doBoostROC, method = method, data = data.original, pauc = pauc, p = p, bw = bw)
            }

        resBoot <- simplify2array(resBoot)    
        rockb <- simplify2array(resBoot["ROC",])
        auckb <- unlist(resBoot["AUC",])
        if(pauc$compute){
            pauckb <- unlist(resBoot["pAUC",])
        }
    }
    columns <-switch(as.character(B > 0), "TRUE" = 1:3, "FALSE" = 1)
    col.names <-c("est","ql", "qh")[columns]
    
    poolROC <- matrix(0, ncol = length(columns), nrow = np, dimnames = list(1:np, col.names))
    poolROC[,1] <- rock
    AUC <- auck
    if(pauc$compute){
        pAUC <- pauck
    }
    if(B > 0) {
        poolROC[,2] <- apply(rockb,1, quantile, prob = 0.025)
        poolROC[,3] <- apply(rockb,1, quantile, prob = 0.975)
        AUC <- c(AUC, quantile(auckb,c(0.025,0.975)))
        if(pauc$compute){
            pAUC <- c(pAUC, quantile(pauckb,c(0.025,0.975)))
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
    res$bws <- obj$bws
    res$bw <- bw
    res$p <- p
    res$ROC <- poolROC
    res$AUC <- AUC
    if(pauc$compute){
        res$pAUC <- if(pauc$focus == "FPF") {
        	pAUC/pauc$value
        } else {
        	pAUC/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    class(res) <- c("pooledROC.kernel", "pooledROC")
    res
}
