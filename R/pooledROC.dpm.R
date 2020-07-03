pooledROC.dpm <-
function(marker, group, tag.h, data,
standardise = TRUE,
p = seq(0,1,l = 101),
compute.lpml = FALSE,
compute.WAIC = FALSE,
compute.DIC = FALSE,
pauc = pauccontrol(),
density = densitycontrol(),
prior.h = priorcontrol.dpm(),
prior.d = priorcontrol.dpm(),
mcmc = mcmccontrol(), 
parallel = c("no", "multicore", "snow"), 
ncpus = 1, 
cl = NULL) {
    
    doMCMCROC <- function(k, res0, res1, L.h, L.d, pauc, density) {
        res <- list()
        if(L.h == 1 & L.d == 1) {
            # Binormal model
            a <- (res0$Mu[k] - res1$Mu[k])/sqrt(res1$Sigma[k])
            b <- sqrt(res0$Sigma[k])/sqrt(res1$Sigma[k])
            
            # ROC curve
            res$ROC <- 1 - pnorm(a + b*qnorm(1-p))
            # AUC
            res$AUC <- 1 - pnorm(a/sqrt(1+b^2))

            if(pauc$compute) {
                if(pauc$focus == "FPF"){
                    res$pAUC <- pbivnorm(-a/sqrt(1+b^2), qnorm(pauc$value), -b/sqrt(1+b^2))
                } else {
                    res$pAUC <- pbivnorm(-a/sqrt(1+b^2), qnorm(1-pauc$value), -1/sqrt(1+b^2))
                }
            }

            if(density$compute) {
                res$dens.h <- dnorm(density$grid.h, mean = res0$Mu[k], sd = sqrt(res0$Sigma[k]))
                res$dens.d <- dnorm(density$grid.d, mean = res1$Mu[k], sd = sqrt(res1$Sigma[k]))
            }
            
        } else if(L.h == 1 & L.d > 1){
            aux1 <- norMix(mu = res1$Mu[k,], sigma = sqrt(res1$Sigma[k,]), w = res1$P[k,])
            q0 <- qnorm(1-p, mean = res0$Mu[k], sd= sqrt(res0$Sigma[k]))
            res$ROC <-  1 - pnorMix(q0, aux1)
            # res$AUC <- simpson(rocdpm[,k], p)
            ##############################################################
            # Computed using results in Erkanly et al (Stat Med, 2006)
            ##############################################################
            a <- outer(res1$Mu[k,], res0$Mu[k], "-")*(sqrt(res1$Sigma[k,])^(-1))
            b <- outer(sqrt(res1$Sigma[k,])^(-1), sqrt(res0$Sigma[k]), "*")
            
            auc_aux <- pnorm(a/sqrt(1+b^2))
            prod.weights <- outer(res1$P[k,], 1, "*")

            res$AUC <- sum(auc_aux*prod.weights)
            ###############################################################
            if(pauc$compute) {
                if(pauc$focus == "FPF"){
                    q01 <- qnorm(1-pauc$pu, mean = res0$Mu[k], sd= sqrt(res0$Sigma[k]))
                    rocdpm1 <-  1 - pnorMix(q01, aux1)
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(pauc$value), -c(b/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                } else {
                    q1 <-qnorMix(1-pauc$pu, aux1)
                    rocdpm1 <-  pnorm(q1, mean = res0$Mu[k], sd = sqrt(res0$Sigma[k]))
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(1-pauc$value), -c(1/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                }
            }

            if(density$compute) {
                res$dens.h <- dnorm(density$grid.h, mean = res0$Mu[k], sd = sqrt(res0$Sigma[k]))
                res$dens.d <- dnorMix(x = density$grid.d, aux1)
            }

        } else if (L.h > 1 & L.d == 1){
            aux0 <- norMix(mu = res0$Mu[k,], sigma = sqrt(res0$Sigma[k,]), w = res0$P[k,])
            q0 <- qnorMix(1-p, aux0)
            res$ROC <-  1 - pnorm(q0, mean = res1$Mu[k], sd = sqrt(res1$Sigma[k]))
            #res$AUC <- simpson(rocdpm[,k], p)
            ##############################################################
            # Computed using results in Erkanly et al (Stat Med, 2006)
            ##############################################################
            a <- outer(res1$Mu[k], res0$Mu[k,], "-")*(sqrt(res1$Sigma[k])^(-1))
            b <- outer(sqrt(res1$Sigma[k])^(-1), sqrt(res0$Sigma[k,]), "*")
            
            auc_aux <- pnorm(a/sqrt(1+b^2))
            prod.weights <- outer(1, res0$P[k,], "*")

            res$AUC <- sum(auc_aux*prod.weights)
            ##############################################################
            if(pauc$compute) {
                if(pauc$focus == "FPF"){
                    q01 <- qnorMix(1-pauc$pu, aux0)
                    rocdpm1 <-  1 - pnorm(q01, mean = res1$Mu[k], sd = sqrt(res1$Sigma[k]))
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(pauc$value), -c(b/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                } else{
                    q1 <- qnorm(1-pauc$pu, mean = res1$Mu[k], sd = sqrt(res1$Sigma[k]))
                    rocdpm1 <-  pnorMix(q1,aux0)
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(1-pauc$value), -c(1/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                }
            }
            if(density$compute) {
                res$dens.h <- dnorMix(x = density$grid.h, aux0)
                res$dens.d <- dnorm(density$grid.d, mean = res1$Mu[k], sd = sqrt(res1$Sigma[k]))                
            }
        } else{
            aux0 <- norMix(mu = res0$Mu[k,], sigma = sqrt(res0$Sigma[k,]), w = res0$P[k,])
            aux1 <- norMix(mu = res1$Mu[k,], sigma = sqrt(res1$Sigma[k,]), w = res1$P[k,])
            q0 <- qnorMix(1-p, aux0)
            res$ROC <- 1 - pnorMix(q0, aux1)
            #res$AUC <- simpson(rocdpm[,k], p)
            ##############################################################
            # Computed using results in Erkanly et al (Stat Med, 2006)
            ###############################################################
            a <- outer(res1$Mu[k,], res0$Mu[k,], "-")*(sqrt(res1$Sigma[k,])^(-1))
            b <- outer(sqrt(res1$Sigma[k,])^(-1), sqrt(res0$Sigma[k,]), "*")
            
            auc_aux <- pnorm(a/sqrt(1+b^2))
            prod.weights <- outer(res1$P[k,], res0$P[k,], "*")

            res$AUC <- sum(auc_aux*prod.weights)
            ###############################################################
            if(pauc$compute) {
                if(pauc$focus == "FPF"){
                    q01 <- qnorMix(1-pauc$pu, aux0)
                    rocdpm1 <- 1 - pnorMix(q01, aux1)
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(pauc$value), -c(b/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                } else {
                    q1 <- qnorMix(1-pauc$pu, aux1)
                    rocdpm1 <- pnorMix(q1, aux0)
                    res$pAUC <- simpson(rocdpm1, pauc$pu)
                    #aux <- pbivnorm(c(a/sqrt(1+b^2)), qnorm(1-pauc$value), -c(1/sqrt(1+b^2)))
                    #res$pAUC <- sum(aux*c(prod.weights))
                }
            }
             if(density$compute) {
                res$dens.h <- dnorMix(x = density$grid.h, aux0)
                res$dens.d <- dnorMix(x = density$grid.d, aux1)                
            }
        }
        res
    }
    parallel <- match.arg(parallel)

    pauc <- do.call("pauccontrol", pauc)    
    density <- do.call("densitycontrol", density)
    mcmc <- do.call("mcmccontrol", mcmc)
    prior.h <- do.call("priorcontrol.dpm", prior.h)
    prior.d <- do.call("priorcontrol.dpm", prior.d)
    
    # Obtain the marker in healthy and diseased
    yh <- data[,marker][data[,group] == tag.h]
    yd <- data[,marker][data[,group] != tag.h]

    # Missing data
    omit.h <- is.na(yh)
    omit.d <- is.na(yd)
    
    yh.wom <- yh[!omit.h]
    yd.wom <- yd[!omit.d]
    
    # Priors
    L.d <- prior.d$L
    m0.d <- prior.d$m0
    S0.d <- prior.d$S0
    a.d <- prior.d$a
    b.d <- prior.d$b
    aalpha.d <- prior.d$aalpha
    balpha.d <- prior.d$balpha
    
    L.h <- prior.h$L
    m0.h <- prior.h$m0
    S0.h <- prior.h$S0
    a.h <- prior.h$a
    b.h <- prior.h$b
    aalpha.h <- prior.h$aalpha
    balpha.h <- prior.h$balpha
    
    if(is.na(L.d)) {
        L.d <- 10 
    } else {
        if(length(L.d) != 1) {
            stop(paste0("L must be a constant"))
        }
    }
    
    if(is.na(L.h)) {
        L.h <- 10 
    } else {
        if(length(L.h) != 1) {
            stop(paste0("L must be a constant"))
        }
    }

    if(is.na(m0.h)) {
        if(standardise == TRUE) m0.h <- 0
        else m0.h <- mean(yh.wom)
    } else { 
        if(length(m0.h) != 1) {
            stop(paste0("m0 must be a constant"))
        }
    }
    
    if(is.na(m0.d)) {
        if(standardise == TRUE) m0.d <- 0
        else m0.d <- mean(yd.wom)
    } else { 
        if(length(m0.d) != 1) {
            stop(paste0("m0 must be a constant"))
        }
    }
    
    if(is.na(S0.h)) {
        if(standardise == TRUE) S0.h <- 10
        else S0.h <- 100
    } else { 
        if(length(S0.h) != 1) {
            stop(paste0("S0 must be a constant"))
        }
    }
    
    if(is.na(S0.d)) {
        if(standardise == TRUE) S0.d <- 10
        else S0.d <- 100
    } else { 
        if(length(S0.d) != 1) {
            stop(paste0("S0 must be a constant"))
        }
    }
    
    if(is.na(a.h)) {
        if(standardise == TRUE) a.h <- 2
        else a.h <- 2
    } else { 
        if(length(a.h) != 1) {
            stop(paste0("a must be a constant"))
        }
    }
    
    if(is.na(b.h)) {
        if(standardise == TRUE) b.h <- 2
        else b.h <- var(yh.wom)
    } else { 
        if(length(b.h) != 1) {
            stop(paste0("b must be a constant"))
        }
    }
    
    if(is.na(a.d)) {
        if(standardise == TRUE) a.d <- 2
        else a.d <- 2
    } else { 
        if(length(a.d) != 1) {
            stop(paste0("a must be a constant"))
        }
    }
    
    if(is.na(b.d)) {
        if(standardise == TRUE) b.d <- 2
        else b.d <- var(yd.wom)
    } else { 
        if(length(b.d) != 1) {
            stop(paste0("b must be a constant"))
        }
    }
    
    if(L.h > 1) {
        if(is.na(aalpha.h)) {
            aalpha.h <- 2 
        } else { 
            if(length(aalpha.h) != 1) {
                stop(paste0("aalpha must be a constant"))
            }
        }
        if(is.na(balpha.h)) {
            balpha.h <- 2 
        } else{ 
            if(length(balpha.h) != 1) {
                stop(paste0("balpha must be a constant"))
            }
        }
    }
    
    if(L.d > 1){
        if(is.na(aalpha.d)) {
            aalpha.d <- 2 
        } else { 
            if(length(aalpha.d) != 1) {
                stop(paste0("aalpha must be a constant"))
            }
        }
        
        if(is.na(balpha.d)) {
            balpha.d <- 2 
        } else { 
            if(length(balpha.d) != 1) {
                stop(paste0("balpha must be a constant"))
            }
        }
    }
    np <- length(p)
    
    # Check if the seq of FPF is correct
    if((L.h > 1 | L.d > 1) & pauc$compute) {
        if(np%%2 == 0) {
            warning("The set of FPFs at which the pooled ROC curve is estimated is not correct. The set is used for calculating the pAUC using Simpson's rule. As such, it should have an odd length.")
        }
    }

    if(L.h > 1){
        res0 <- dpm(y = yh.wom,
        prior = list(m0 = m0.h,
        S0 = S0.h,
        a = a.h,
        b = b.h,
        aalpha = aalpha.h,
        balpha = balpha.h,
        L = L.h),
        mcmc = mcmc,
        standardise = standardise)
    }
    if(L.h == 1){
        res0 <- gibbsn(y = yh.wom,
        prior = list(m0 = m0.h,
        S0 = S0.h,
        a = a.h,
        b = b.h),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    if(L.d > 1){
        res1 <- dpm(y = yd.wom,
        prior = list(m0 = m0.d,
        S0 = S0.d,
        a = a.d,
        b = b.d,
        aalpha = aalpha.d,
        balpha = balpha.d,
        L = L.d),
        mcmc = mcmc,
        standardise = standardise)
    }
    if(L.d == 1){
        res1 <- gibbsn(y = yd.wom,
        prior = list(m0 = m0.d,
        S0 = S0.d,
        a = a.d,
        b = b.d),
        mcmc = mcmc,
        standardise = standardise)
    }

    if(density$compute){ 
        if(all(is.na(density$grid.h))) {
            density$grid.h <- seq(min(yh.wom) - 1, max(yh.wom) + 1, len = 200)
        }        
        if(all(is.na(density$grid.d))) {
            density$grid.d <- seq(min(yd.wom) - 1, max(yd.wom) + 1, len = 200)
        }
    }
    if(pauc$compute) {
        if(pauc$focus == "FPF") {
            pauc$pu <- seq(0, pauc$value, len = np)
        } else {
            pauc$pu <- seq(pauc$value, 1, len = np)
        }
    }    
    if(mcmc$nsave > 0) {
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
                    parallel::mclapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, pauc = pauc, density = density, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, pauc = pauc, density = density)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, pauc = pauc, density = density)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, pauc = pauc, density = density)
            }

        resBoot <- simplify2array(resBoot)    
        rocdpm <- simplify2array(resBoot["ROC",])
        aucdpm <- unlist(resBoot["AUC",])
        if(pauc$compute){
            paucdpm <- unlist(resBoot["pAUC",])
        }
        if(density$compute){
            dens.h <- t(simplify2array(resBoot["dens.h",]))
            dens.d <- t(simplify2array(resBoot["dens.d",]))
        }

    } else {
        stop("nsave should be larger than zero.")
    }
    
    poolROC <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
    poolROC[,1] <- apply(rocdpm, 1, mean)
    poolROC[,2] <- apply(rocdpm, 1, quantile, prob = 0.025)
    poolROC[,3] <- apply(rocdpm, 1, quantile, prob = 0.975)
    
    res <- list()
    res$call <- match.call()
    res$marker <- list(h = yh, d = yd)
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$p <- p
    res$mcmc <- mcmc
    res$prior <- list()
    if(L.d == 1){
        res$prior$d <- list(m0 = m0.d,
        S0 = S0.d,
        a = a.d,
        b = b.d,
        L = L.d)
    }
    if(L.d > 1){
        res$prior$d <- list(m0 = m0.d,
        S0 = S0.d,
        a = a.d,
        b = b.d,
        aalpha = aalpha.d,
        balpha = balpha.d,
        L = L.d)
    }
    if(L.h == 1){
        res$prior$h <- list(m0 = m0.h,
        S0 = S0.h,
        a = a.h,
        b = b.h,
        L = L.h)
    }
    if(L.h > 1){
        res$prior$h <- list(m0 = m0.h,
        S0 = S0.h,
        a = a.h,
        b = b.h,
        aalpha = aalpha.h,
        balpha = balpha.h,
        L = L.h)
    }
    res$ROC <- poolROC
    AUC <- c(mean(aucdpm), quantile(aucdpm, c(0.025,0.975)))
    names(AUC) <- c("est","ql", "qh")
    res$AUC <- AUC
    
    if(pauc$compute) {
        pAUC <- c(mean(paucdpm), quantile(paucdpm, c(0.025,0.975)))
        names(pAUC) <- c("est","ql", "qh")
        res$pAUC <- if(pauc$focus == "FPF") {
            pAUC/pauc$value
        } else {
            pAUC/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    if(density$compute){
        res$dens <- list()
        res$dens$h <- list(grid = density$grid.h, dens = dens.h)
        res$dens$d <- list(grid = density$grid.d, dens = dens.d)
    }
    if(compute.lpml | compute.WAIC | compute.DIC) {
        termsumh <- inf_criteria_dpm(y = yh.wom, res = res0, L = L.h)
        termsumd <- inf_criteria_dpm(y = yd.wom, res = res1, L = L.d)
    }
    if(compute.lpml) {
        res$lpml <- list()
        res$lpml$h <- lpml_dpm(y = yh.wom, res = res0, L = L.h, termsum = termsumh)
        res$lpml$d <- lpml_dpm(y = yd.wom, res = res1, L = L.d, termsum = termsumd)
    }
    if(compute.WAIC) {
        res$WAIC <- list()
        res$WAIC$h <- waic_dpm(y = yh.wom, res = res0, L = L.h, termsum = termsumh)
        res$WAIC$d <- waic_dpm(y = yd.wom, res = res1, L = L.d, termsum = termsumd)
    }
    if(compute.DIC) {
        res$DIC <- list()
        res$DIC$h <- dic_dpm(y = yh.wom, res = res0, L = L.h, termsum = termsumh)
        res$DIC$d <- dic_dpm(y = yd.wom, res = res1, L = L.d, termsum = termsumd)
    }
    
    #res$fit <- list(h = list(P = p0, Mu = mu0, Sigma2 = sigma02), d = list(P = p1, Mu = mu1, Sigma2 = sigma12))
    res$fit <- list(h = res0, d = res1)
    class(res) <- c("pooledROC.dpm", "pooledROC")
    res
}
