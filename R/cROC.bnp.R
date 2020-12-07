cROC.bnp <-
function(formula.h,
formula.d,
group,
tag.h,
data,
newdata,
standardise = TRUE,
p = seq(0, 1, l = 101),
compute.lpml = FALSE,
compute.WAIC = FALSE,
compute.DIC = FALSE,
pauc = pauccontrol(),
density = densitycontrol(),
prior.h = priorcontrol.bnp(),
prior.d = priorcontrol.bnp(),
mcmc = mcmccontrol(), 
parallel = c("no", "multicore", "snow"), 
ncpus = 1, 
cl = NULL) {
    doMCMCROC <- function(k, res0, res1, L.h, L.d, Xhp, Xdp, pauc, density, p) {
        np <- length(p)
        npred <- nrow(Xhp)

        npred <- nrow(newdata)
        rocddp <- matrix(0, nrow = np, ncol = npred)
        aucddp <- vector(length = npred)
        
        meanfun.h <- meanfun.d <- vector(length = npred)

        if(pauc$compute) {
            rocddp_pauc <- matrix(0, nrow = np, ncol = npred)
            paucddp <- vector(length = npred)
        }

        if(density$compute){
            dens.h <- matrix(0, nrow = length(density$grid.h), ncol = npred)
            dens.d <- matrix(0, nrow = length(density$grid.d), ncol = npred)
        }


        if(L.d == 1 & L.h == 1){
            mu.h <- Xhp%*%res0$Beta[k,]
            mu.d <- Xdp%*%res1$Beta[k,]
            
            meanfun.d <- mu.d
            meanfun.h <- mu.h
            
            # Binormal model
            a <- as.vector(mu.h - mu.d)/sqrt(res1$Sigma2[k])
            b <- sqrt(res0$Sigma2[k])/sqrt(res1$Sigma2[k]) # ROC curve
            rocddp <- t(1 - apply(a + outer(rep(b, npred), qnorm(1-p), "*"), c(1,2), pnorm))
            
            # AUC
            aucddp <- 1 - pnorm(a/sqrt(1+b^2))

            if(pauc$compute) {                    
                if(pauc$focus == "FPF") {
                    paucddp <- pbivnorm(-a/sqrt(1+b^2), qnorm(pauc$value), -b/sqrt(1+b^2))
                } else {
                    paucddp <- pbivnorm(-a/sqrt(1+b^2), qnorm(1-pauc$value), -1/sqrt(1+b^2))
                }
            }

            for(l in 1:npred){                
                if(density$compute){
                    dens.h[,l] <- dnorm(density$grid.h, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                    dens.d[,l] <- dnorm(density$grid.d, mean = mu.d[l], sd = sqrt(res1$Sigma2[k]))
                }
            }
        }        
        if(L.d == 1 & L.h > 1) {
            mu.h <- tcrossprod(Xhp, res0$Beta[k,,])
            mu.d <- Xdp%*%res1$Beta[k,]
            meanfun.d <- mu.d
            
            for(l in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = sqrt(res0$Sigma2[k,]), w = res0$P[k,])
                q0 <- qnorMix(1 - p, aux0)
                rocddp[,l] <- 1 - pnorm(q0, mean = mu.d[l], sd = sqrt(res1$Sigma2[k]))
                aucddp[l] <- simpson(rocddp[,l], p)
                meanfun.h[l] = sum(res0$P[k,]*t(mu.h[l,]))
                
                if(pauc$compute) {      
                    if(pauc$focus == "FPF"){
                        pu <- seq(p[1], pauc$value, len = np)
                        q0_pauc <- qnorMix(1-pu, aux0)
                        rocddp_pauc[,l] <- 1 - pnorm(q0_pauc, mean= mu.d[l,], sd = sqrt(res1$Sigma2[k]))
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    } else{
                        pu <- seq(pauc$value, p[np], len = np)
                        q1_pauc <- qnorm(1-pu, mean= mu.d[l], sd = sqrt(res1$Sigma2[k]))
                        rocddp_pauc[,l] <- pnorMix(q1_pauc, aux0)
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    }
                }
                
                if(density$compute){
                    dens.d[,l] <- dnorm(density$grid.d, mean = mu.d[l], sd = sqrt(res1$Sigma2[k]))
                    dens.h[,l] <- dnorMix(density$grid.h, aux0)
                }
            }
        }
        
        if(L.d > 1 & L.h == 1) {
            mu.h <- Xhp%*%res0$Beta[k,]
            mu.d <- tcrossprod(Xdp, res1$Beta[k,,]) #Xdp%*%t(res1$Beta[k,,])
            meanfun.h <- mu.h
            
            for(l in 1:npred) {
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = sqrt(res1$Sigma2[k,]), w = res1$P[k,])
                q0 <- qnorm(1-p, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                rocddp[,l] <- 1 - pnorMix(q0, aux1)
                aucddp[l] <- simpson(rocddp[,l], p)
                
                meanfun.d[l] = sum(res0$P[k,]*t(mu.d[l,]))
                
                if(pauc$compute) {
                    
                    if(pauc$focus == "FPF"){
                        pu <- seq(p[1], pauc$value, len = np)
                        q0_pauc <- qnorm(1-pu, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                        rocddp_pauc[,l] <- 1 - pnorMix(q0_pauc, aux1)
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    } else{
                        pu <- seq(pauc$value, p[np], len = np)
                        q1_pauc <- qnorMix(1-pu, aux1)
                        rocddp_pauc[,l] <- pnorm(q1_pauc, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    }
                }
                if(density$compute){
                    dens.h[,l] <- dnorm(density$grid.h, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                    dens.d[,l] <- dnorMix(density$grid.d, aux1)
                }
            }
        }
        
        if(L.d > 1 & L.h > 1) {
            mu.h <- tcrossprod(Xhp, res0$Beta[k,,])
            mu.d <- tcrossprod(Xdp, res1$Beta[k,,])
            
            for(l in 1:npred) {
                aux0 <- norMix(mu = c(mu.h[l,]), sigma = sqrt(res0$Sigma2[k,]), w = res0$P[k,])
                aux1 <- norMix(mu = c(mu.d[l,]), sigma = sqrt(res1$Sigma2[k,]), w = res1$P[k,])
                
                q0 <- qnorMix(1-p, aux0)
                rocddp[,l] <- 1 - pnorMix(q0, aux1)
                aucddp[l] <- simpson(rocddp[,l], p)
                
                meanfun.h[l] <- sum(res0$P[k,]*t(mu.h[l,]))
                meanfun.d[l] <- sum(res1$P[k,]*t(mu.d[l,]))
                
                if(pauc$compute) { 
                    if(pauc$focus == "FPF"){
                        pu <- seq(p[1], pauc$value, len = np)
                        q0_pauc <- qnorMix(1-pu, aux0)
                        rocddp_pauc[,l] <- 1 - pnorMix(q0_pauc, aux1)
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    } else{
                        pu <- seq(pauc$value, p[np], len = np)
                        q1_pauc <- qnorMix(1-pu, aux1)
                        rocddp_pauc[,l] <- pnorMix(q1_pauc, aux0)
                        paucddp[l] <- simpson(rocddp_pauc[,l], pu)
                    }
                }
                if(density$compute){
                    dens.d[,l] <- dnorMix(density$grid.d, aux1)
                    dens.h[,l] <- dnorMix(density$grid.h, aux0)
                }
            }
        }

        res <- list()
        res$ROC <- rocddp
        res$AUC <- aucddp
        res$meanfun.h <- meanfun.h
        res$meanfun.d <- meanfun.d
        if(pauc$compute) {
            res$pAUC <- paucddp
        }
        if(density$compute) {
            res$dens.h <- dens.h
            res$dens.d <- dens.d
        }
        res
    }
    parallel <- match.arg(parallel)

    pauc <- do.call("pauccontrol", pauc)
    density <- do.call("densitycontrol", density)
    mcmc <- do.call("mcmccontrol", mcmc)
    prior.h <- do.call("priorcontrol.bnp", prior.h)
    prior.d <- do.call("priorcontrol.bnp", prior.d)

    if(inherits(formula.h, "character")) {
        formula.h <- as.formula(formula.h)
    }
    if(inherits(formula.d, "character")) {
        formula.d <- as.formula(formula.d)
    }
    
    # Marker variable
    tf <- terms.formula(formula.h, specials = c("f"))
    if (attr(tf, "response") > 0) {
        marker.h <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The 'formula.h' should include the response variable (left hand side)")
    }
    
    tf <- terms.formula(formula.d, specials = c("f"))
    if (attr(tf, "response") > 0) {
        marker.d <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The 'formula.h' should include the response variable (left hand side)")
    }
    
    if(marker.h != marker.d) {
        stop("The response variable (biomarker) in 'formula.h' and 'formula.d' should be the same")
    }

    marker <- marker.h
    
    # Variables in the model
    names.cov.h <- get_vars_formula(formula.h) #all.vars(formula.h)[-1]
    names.cov.d <- get_vars_formula(formula.d) #all.vars(formula.d)[-1]
    names.cov <- c(names.cov.h, names.cov.d[is.na(match(names.cov.d, names.cov.h))])
    
    if(!missing(newdata) && !inherits(newdata, "data.frame"))
        stop("Newdata must be a data frame")
    if(sum(is.na(match(c(marker,names.cov,group), names(data)))))
        stop("Not all needed variables are supplied in data")
    if(!missing(newdata) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(newdata)))))
        stop("Not all needed variables are supplied in newdata")
    if(length(unique(data[,group]))!= 2)
        stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep="")
     
    # New data, removing missing values
    data.new <- data[,c(marker,group,names.cov)]
    omit.h <- apply(data.new[data.new[,group] == tag.h, c(marker, group, names.cov.h)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.h, c(marker, group, names.cov.d)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.h,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.h,,drop = FALSE][!omit.d,,drop = FALSE])
    
    data.h <- data.new[data.new[,group] == tag.h,]
    data.d <- data.new[data.new[,group] != tag.h,]
    
    n0 <- nrow(data.h)
    n1 <- nrow(data.d)
    np <- length(p)
    
    # Construct design matrix: healthy
    MM0 <- design.matrix.bnp(formula.h, data.h, standardise)
    X0 <- MM0$X
    k0 <- ncol(X0)
    
    # Construct design matrix: diseased
    MM1 <- design.matrix.bnp(formula.d, data.d, standardise)
    X1 <- MM1$X
    k1 <- ncol(X1)
    
    # New data (for predictions)
    if(missing(newdata)) {
        newdata <- cROCData(data.new, names.cov, group)
    } else {
        newdata <- na.omit(newdata[,names.cov,drop=FALSE])
    }
    
    data.h.marker <- data.h[,marker]
    data.d.marker <- data.d[,marker]

    # Getting OLS estimates
    #coefs.h <- solve(t(X0) %*% X0) %*% t(X0) %*% data.h.marker
    #var.h <- sum((data.h.marker - X0 %*% coefs.h)^2)/(n0 - ncol(X0))
    #cov.h <- solve(t(X0) %*% X0)*var.h
    res <- ols.function(X0, data.h.marker, vcov = TRUE)
    coefs.h <- res$coeff
    var.h <- sum((data.h.marker - X0 %*% coefs.h)^2)/(n0 - ncol(X0))
    cov.h <- res$vcov*var.h
    
    #coefs.d <- solve(t(X1) %*% X1) %*% t(X1) %*% data.d.marker
    #var.d <- sum((data.d.marker - X1 %*% coefs.d)^2)/(n1 - ncol(X1))
    #cov.d <- solve(t(X1) %*% X1)*var.d
    res <- ols.function(X1, data.d.marker, vcov = TRUE)
    coefs.d <- res$coeff
    var.d <- sum((data.d.marker - X1 %*% coefs.d)^2)/(n1 - ncol(X1))
    cov.d <- res$vcov*var.d
    
    # Hyperparameters
    L.d <- prior.d$L
    m0.d <- prior.d$m0
    S0.d <- prior.d$S0
    nu.d <- prior.d$nu
    Psi.d <- prior.d$Psi
    a.d <- prior.d$a
    b.d <- prior.d$b
    aalpha.d <- prior.d$aalpha
    balpha.d <- prior.d$balpha
    
    L.h <- prior.h$L
    m0.h <- prior.h$m0
    S0.h <- prior.h$S0
    nu.h <- prior.h$nu
    Psi.h <- prior.h$Psi
    a.h <- prior.h$a
    b.h <- prior.h$b
    aalpha.h <- prior.h$aalpha
    balpha.h <- prior.h$balpha
    
    if(is.na(L.h)) {
        L.h <- 10 
    } else { 
        if(length(L.h) != 1) {
            stop(paste0("L must be a constant"))
        }
    }
    
    if(is.na(L.d)) {
        L.d <- 10 
    } else {
        if(length(L.d) != 1) {
            stop(paste0("L must be a constant"))
        }
    }
        
    if(is.na(m0.h)) {
        if(standardise) m0.h <- rep(0, k0)
        else m0.h <- coefs.h
    } else { 
        if(length(m0.h) != k0) {
            stop(paste0("'m0.h' must be a vector of length ", k0))
        }
    }
    
    if(is.na(m0.d)) {
        if(standardise) m0.d <- rep(0, k1)
        else m0.d <- coefs.d
    } else { 
        if(length(m0.d) != k1) {
            stop(paste0("'m0.d' must be a vector of length ", k1))
        }
    }
    
    if(is.na(S0.h)) {
        if(standardise) S0.h <- 10*diag(k0)
        else S0.h <- cov.h
    } else { 
        if(!is.matrix(S0.h) | !all(dim(S0.h) == c(k0, k0))) {
            stop(paste0("'S0.h' must be a matrix of dimension ", k0, "x", k0))
        }
    }
    
    if(is.na(S0.d)) {
        if(standardise) S0.d <- 10*diag(k1)
        else S0.d <- cov.d
    } else { 
        if(!is.matrix(S0.d) | !all(dim(S0.d) == c(k1, k1))) {
            stop(paste0("'S0.d' must be a matrix of dimension ", k1, "x", k1))
        }
    }
    
    if(is.na(Psi.h)) {
        if(standardise) Psi.h <- diag(k0)
        else Psi.h <- 30*cov.h
    } else { 
        if(!is.matrix(Psi.h) | !all(dim(Psi.h) == c(k0, k0))) {
            stop(paste0("'Psi.h' must be a matrix of dimension ", k0, "x", k0))
        }
    }
    
    if(is.na(Psi.d)) {
        if(standardise) Psi.d <- diag(k1)
        else Psi.d <- 30*cov.d
    } else { 
        if(!is.matrix(Psi.d) | !all(dim(Psi.d) == c(k1, k1))) {
            stop(paste0("'Psi.d' must be a matrix of dimension ", k1, "x", k1))
        }
    }
    
    if(is.na(nu.h)) {
        nu.h <- k0 + 2
    } else { 
        if(nu.h < k0 + 2){
            stop(paste0("'nu.h' must be larger than ", k0 + 2))
        }
    }
    
    if(is.na(nu.d)) { 
        nu.d <- k1 + 2
    } else { 
        if(nu.d < k1 + 2){
            stop(paste0("'nu.h' must be larger than ", k1 + 2))
        }
    }
    
    if(is.na(a.h)) {
        a.h <- 2
    } else { 
        if(length(a.h) != 1) {
            stop(paste0("'a.h' must be a constant"))
        }
    }
    
    if(is.na(a.d)) { 
        a.d <- 2
    } else { 
        if(length(a.d) != 1) {
            stop(paste0("'a.d' must be a constant"))
        }
    }
    
    if(is.na(b.h)) {
        if(standardise) b.h <- 2
        else b.h <- var.h
    } else { 
        if(length(b.h) != 1) {
            stop(paste0("'b.h' must be a constant"))
        }
    }
    
    if(is.na(b.d)) {
        if(standardise) b.d <- 2
        else b.d <- var.d
    } else { 
        if(length(b.d) != 1) {
            stop(paste0("'b.d' must be a constant"))
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
        } else { 
            if(length(balpha.h) != 1) {
                stop(paste0("balpha must be a constant"))
            }
        }
    }
    
    if(L.d > 1) {
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
    
    if(L.h > 1 | L.d > 1) {
        if(min(p) != 0 | max(p) != 1 | np%%2 == 0) {
            warning("The set of FPFs at which the covariate-specific ROC curve is estimated is not correct. The set is used for calculating the AUC using Simpson's rule. As such, it should (a) start in 0; (b) end in 1; and (c) have an odd length.")
        }
    }

    if(L.h == 1){
        res0 <- regnth(y = data.h.marker,
        X = X0,
        prior = list(m0 = m0.h,
        S0 = S0.h,
        nu = nu.h,
        Psi = Psi.h,
        a = a.h,
        b = b.h),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    if(L.h > 1){
        res0 <- bddp(y = data.h.marker,
        X = X0,
        prior = list(m0 = m0.h,
        S0 = S0.h,
        nu = nu.h,
        Psi = Psi.h,
        a = a.h,
        b = b.h,
        aalpha = aalpha.h,
        balpha = balpha.h,
        L = L.h),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    if(L.d == 1){
        res1 <- regnth(y = data.d.marker, X = X1,
        prior = list(m0 = m0.d, S0 = S0.d,
        nu = nu.d, Psi = Psi.d,
        a = a.d, b = b.d),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    if(L.d > 1){
        res1 <- bddp(y = data.d.marker, X = X1,
        prior = list(m0 = m0.d, S0 = S0.d,
        nu = nu.d, Psi = Psi.d,
        a = a.d, b = b.d,
        aalpha = aalpha.d, balpha = balpha.d,
        L = L.d),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    # Compute conditional ROC and AUC
    X0p <- predict(MM0, newdata = newdata)$X
    X1p <- predict(MM1, newdata = newdata)$X
    
    if(density$compute){
        if(all(is.na(density$grid.h))) {
            density$grid.h <- seq(min(data.h.marker) - 1, max(data.h.marker) + 1, len = 200)
        }
        
        if(all(is.na(density$grid.d))) {
            density$grid.d <- seq(min(data.d.marker) - 1, max(data.d.marker) + 1, len = 200)
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
                    parallel::mclapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, Xhp = X0p, Xdp = X1p, pauc = pauc, density = density, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, Xhp = X0p, Xdp = X1p, pauc = pauc, density = density, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, Xhp = X0p, Xdp = X1p, pauc = pauc, density = density, p = p)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, res1 = res1, L.h = L.h, L.d = L.d, Xhp = X0p, Xdp = X1p, pauc = pauc, density = density, p = p)
            }

        resBoot <- simplify2array(resBoot)    
        rocddp <- aperm(simplify2array(resBoot["ROC",]), c(1,3,2))
        aucddp <- t(simplify2array(resBoot["AUC",]))
        if(pauc$compute){
            paucddp <- t(simplify2array(resBoot["pAUC",]))
        }
        meanfun.h <- simplify2array(resBoot["meanfun.h",])
        meanfun.d <- simplify2array(resBoot["meanfun.d",])
        if(density$compute){
            dens.h <- aperm(simplify2array(resBoot["dens.h",]), c(1,3,2))
            dens.d <- aperm(simplify2array(resBoot["dens.d",]), c(1,3,2))
        }

    } else {
        stop("nsave should be larger than zero.")
    }

    npred <- nrow(X0p)

    aucddpm <- apply(aucddp, 2, mean)
    aucddpl <- apply(aucddp, 2, quantile, prob = 0.025)
    aucddph <- apply(aucddp, 2, quantile, prob = 0.975)
    
    rocddpm <- rocddpl <- rocddph <- matrix(0, nrow = np, ncol = npred)
    for(l in 1:npred){
        for(i in 1:np){
            rocddpm[i,l] <- mean(rocddp[i,,l])
            rocddpl[i,l] <- quantile(rocddp[i,,l],0.025)
            rocddph[i,l] <- quantile(rocddp[i,,l],0.975)
        }
    }
    
    if(pauc$compute) {
        paucddpm <- apply(paucddp, 2, mean)
        paucddpl <- apply(paucddp, 2, quantile, prob = 0.025)
        paucddph <- apply(paucddp, 2, quantile, prob = 0.975)
    }
    
    meanfun.d.m <- apply(meanfun.d, 1, mean)
    meanfun.d.l <- apply(meanfun.d, 1, quantile, prob = 0.025)
    meanfun.d.h <- apply(meanfun.d, 1, quantile, prob = 0.975)
    
    meanfun.h.m <- apply(meanfun.h, 1, mean)
    meanfun.h.l <- apply(meanfun.h, 1, quantile, prob = 0.025)
    meanfun.h.h <- apply(meanfun.h, 1, quantile, prob = 0.975)
    
    res <- list()
    res$call <- match.call()
    res$newdata <- newdata
    res$data <- data
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$marker <- marker
    res$group <- group
    res$tag.h <- tag.h
    res$mcmc <- mcmc
    res$p <- p
    res$prior <- list()
    if(L.d > 1){
        res$prior$d <-list(m0 = m0.d, S0 = S0.d,
        nu = nu.d, Psi = Psi.d,
        a = a.d, b = b.d,
        aalpha = aalpha.d, balpha = balpha.d,
        L = L.d)
    }
    if(L.d == 1){
        res$prior$d <-list(m0 = m0.d, S0 = S0.d,
        nu = nu.d, Psi = Psi.d,
        a = a.d, b = b.d,
        L = L.d)
    }
    if(L.h > 1){
        res$prior$h <- list(m0 = m0.h, S0 = S0.h,
        nu = nu.h, Psi = Psi.h,
        a = a.h, b = b.h,
        aalpha = aalpha.h, balpha = balpha.h,
        L = L.h)
    }
    if(L.h == 1){
        res$prior$h <- list(m0 = m0.h, S0 = S0.h,
        nu = nu.h, Psi = Psi.h,
        a = a.h, b = b.h,
        L = L.h)
    }
    res$ROC <- list(est = t(rocddpm), ql = t(rocddpl), qh = t(rocddph))
    res$AUC <- data.frame(AUC = aucddpm, AUCql = aucddpl, AUCqh = aucddph)
    
    if(pauc$compute) {
        #res$pAUC <- data.frame(pAUC = paucddpm/pauc$value, pAUCql = paucddpl/pauc$value, pAUCqh = paucddph/pauc$value)
        aux <- data.frame(pAUC = paucddpm, pAUCql = paucddpl, pAUCqh = paucddph)
        res$pAUC <- if(pauc$focus == "FPF") {
            aux/pauc$value
        } else {
            aux/(1-pauc$value)
        }
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    if(density$compute) {
        res$dens <- list()
        res$dens$h <- list(grid = density$grid.h, dens = dens.h)
        res$dens$d <- list(grid = density$grid.d, dens = dens.d)
    }
    
    res$reg.fun <- list()
    res$reg.fun$h <- data.frame(est = meanfun.h.m, ql = meanfun.h.l, qh = meanfun.h.h)
    res$reg.fun$d <- data.frame(est = meanfun.d.m, ql = meanfun.d.l, qh = meanfun.d.h)
    
    if(compute.lpml | compute.WAIC | compute.DIC) {
        termh <- inf_criteria(y = data.h.marker, X = X0, res = res0)        
        termd <- inf_criteria(y = data.d.marker, X = X1, res = res1)
    }
    
    if(compute.lpml) {
        res$lpml <- list()
        if(L.h > 1){
            res$lpml$h <- lpml(y = data.h.marker, X = X0, res = res0,
            L = L.h, termsum = termh)
        }
        
        if(L.h == 1){
            res$lpml$h <- lpmlp(y = data.h.marker, X = X0, res = res0, term = termh)
        }
        
        if(L.d > 1){
            res$lpml$d <- lpml(y = data.d.marker, X = X1, res = res1,
            L = L.d, termsum = termd)
        }
        
        if(L.d == 1){
            res$lpml$d <- lpmlp(y = data.d.marker, X = X1, res = res1, term = termd)
        }
    }
    if(compute.WAIC) {
        res$WAIC <- list()
        if(L.h > 1) {
            res$WAIC$h <- waicnp(y = data.h.marker, X = X0, res = res0,
            L = L.h, termsum = termh)
        }
        if(L.h == 1) {
            res$WAIC$h <- waicp(y = data.h.marker, X = X0, res = res0, term = termh)
        }
        
        if(L.d > 1) {
            res$WAIC$d <- waicnp(y = data.d.marker, X = X1, res = res1,
            L = L.d, termsum = termd)
        }
        if(L.d == 1) {
            res$WAIC$d <- waicp(y = data.d.marker, X = X1, res = res1, term = termd)
        }
    }
    if(compute.DIC) {
        res$DIC <- list()
        if(L.h > 1){
            res$DIC$h <- dic(y = data.h.marker, X = X0, res = res0,
            L = L.h, termsum = termh)
        }
        if(L.h == 1) {
            res$DIC$h <- dicp(y = data.h.marker, X = X0, res = res0, term = termh)
        }
        if(L.d > 1) {
            res$DIC$d <- dic(y = data.d.marker, X = X1, res = res1,
            L = L.d, termsum = termd)
        }
        if(L.d == 1) {
            res$DIC$d <- dicp(y = data.d.marker, X = X1, res = res1, term = termd)
        }
    }
    res$fit <- list()
    if(L.h >1){
        res$fit$h <- list(formula = formula.h,
        mm = MM0,
        beta = res0$Beta,
        sd   = sqrt(res0$Sigma2),
        probs = res0$P)
    }
    
    if(L.h == 1){
        res$fit$h <- list(formula = formula.h,
        mm = MM0,
        beta = res0$Beta,
        sd   = sqrt(res0$Sigma2))
    }
    
    if(L.d > 1){
        res$fit$d <- list(formula = formula.d,
        mm = MM1,
        beta = res1$Beta,
        sd   = sqrt(res1$Sigma2),
        probs = res1$P)
    }
    
    if(L.d == 1){
        res$fit$d <- list(formula = formula.d,
        mm = MM1,
        beta = res1$Beta,
        sd   = sqrt(res1$Sigma2))
    }
    
    res$data_model <- list(y = list(h = data.h.marker,
    d = data.d.marker),
    X = list(h = X0, d = X1))
    res$ci.fit <- TRUE
    class(res) <- c("cROC.bnp", "cROC")
    res
}
