AROC.bnp <-
function(formula.h,
group,
tag.h,
data,
standardise = TRUE,
p = seq(0,1,l = 101),
compute.lpml = FALSE,
compute.WAIC = FALSE,
compute.DIC = FALSE,
pauc = pauccontrol(),
density = densitycontrol.aroc(),
prior.h = priorcontrol.bnp(),
mcmc = mcmccontrol(), 
parallel = c("no", "multicore", "snow"), 
ncpus = 1, 
cl = NULL) {
    
    doMCMCROC <- function(k, res0, L, yd, Xd, pauc, density, p) {
        nd <- length(yd)
        np <- length(p)

        # Placement values
        if(L == 1) {
            up <- 1 - pnorm(yd, mean = Xd%*%res0$Beta[k,], sd = sqrt(res0$Sigma2[k]))
            if(density$compute) {
                # Densities (if needed)
                npred <- nrow(density$Xhp)
                dens.h <- matrix(0, nrow = length(density$grid.h), ncol = npred)
                meanfun.h <- vector(length = npred)

                mu.h <- density$Xhp%*%res0$Beta[k,]
                meanfun.h <- mu.h
                for(l in 1:npred){
                    dens.h[,l] <- dnorm(density$grid.h, mean = mu.h[l], sd = sqrt(res0$Sigma2[k]))
                }
            }
        } else {
            up <- 1 - apply(t(res0$P[k,]*t(pnorm(yd, mean = tcrossprod(Xd, res0$Beta[k,,]), sd = rep(sqrt(res0$Sigma2[k,]), each = length(yd))))), 1, sum)
            if(density$compute) {
                # Densities (if needed)
                npred <- nrow(density$Xhp)
                dens.h <- matrix(0, nrow = length(density$grid.h), ncol = npred)
                meanfun.h <- vector(length = npred)

                mu.h <- tcrossprod(density$Xhp, res0$Beta[k,,])
                for(l in 1:npred){
                    aux0 <- norMix(mu = c(mu.h[l,]), sigma = sqrt(res0$Sigma2[k,]), w = res0$P[k,])
                    dens.h[,l] <- dnorMix(density$grid.h, aux0)
                    meanfun.h[l] <- sum(res0$P[k,]*t(mu.h[l,]))
                }
            }
        }

        aux <- rexp(nd, 1)
        weights <- aux/sum(aux)

        arocddp <- apply(weights*outer(up, p, "<="), 2, sum)
        aaucddp <- 1 - sum(weights*up)

        if(pauc$compute) {
            if(pauc$focus == "FPF"){
                paucddp <- pauc$value - sum(weights*pmin(pauc$value, up))
            } else{
                arocp[1] <- 0
                arocp[np] <- 1
                rocapp <- approxfun(p, arocp, method = "linear")
                p1 <- uniroot(function(x) {rocapp(x) - pauc$value}, interval = c(0, 1))$root
                paucddp <- integrate(rocapp, lower = p1, upper = 1,
                stop.on.error = FALSE)$value - (1 - p1)*pauc$value
            }
        }

        res <- list()
        res$ROC <- arocddp
        res$AUC <- aaucddp
        if(pauc$compute) {
            res$pAUC <- paucddp
        }
        if(density$compute) {
            res$dens.h <- dens.h
            res$meanfun.h <- meanfun.h
        }
        res
    }

    parallel <- match.arg(parallel)

    pauc <- do.call("pauccontrol", pauc)    
    density <- do.call("densitycontrol.aroc", density)
    mcmc <- do.call("mcmccontrol", mcmc)
    prior.h <- do.call("priorcontrol.bnp", prior.h)
    
    if(inherits(formula.h, "character")) {
        formula.h <- as.formula(formula.h)
    }
    # Marker variable
    tf <- terms.formula(formula.h, specials = c("f"))
    if (attr(tf, "response") > 0) {
        marker <- as.character(attr(tf, "variables")[2])
    } else {
        stop("The formula should include the response variable (left hand side)")
    }
    
    # Variables in the model
    names.cov <- get_vars_formula(formula.h) #all.vars(formula.h)[-1]
    
    if(sum(is.na(match(c(marker, names.cov, group), names(data)))))
        stop("Not all needed variables are supplied in data")
    if(length(unique(data[,group])) != 2)
        stop(paste(group," variable must have only two different values (for healthy and diseased individuals)"), sep = "")
    
    # New data, removing missing values
    data.new <- data[,c(marker,group,names.cov)]
    omit.h <- apply(data.new[data.new[,group] == tag.h, c(marker, group, names.cov)], 1, anyNA)
    omit.d <- apply(data.new[data.new[,group] != tag.h, c(marker, group, names.cov)], 1, anyNA)
    
    data.new <- rbind(data.new[data.new[,group] == tag.h,,drop = FALSE][!omit.h,,drop = FALSE], data.new[data.new[,group] != tag.h,,drop = FALSE][!omit.d,,drop = FALSE])
        
    data.h <- data.new[data.new[,group] == tag.h,]
    data.d <- data.new[data.new[,group] != tag.h,]
    
    n0 <- nrow(data.h)
    n1 <- nrow(data.d)
    np <- length(p)
    
    # Construct design matrix
    MM0 <- design.matrix.bnp(formula.h, data.h, standardise)
    X0 <- MM0$X
    k <-  ncol(X0)
    
    # Construct design matrix in diseased population (based on healthy)
    X1 <- predict(MM0, data.d)$X
    
    data.h.marker <- data.h[,marker]
    data.d.marker <- data.d[,marker]

    # Getting OLS estimates
    res <- ols.function(X0, data.h.marker, vcov = TRUE)
    coefs.h <- res$coeff
    var.h <- sum((data.h.marker - X0 %*% coefs.h)^2)/(n0 - ncol(X0))
    cov.h <- res$vcov*var.h
    
    # Hyperparameters
    L <- prior.h$L
    m0 <- prior.h$m0
    S0 <- prior.h$S0
    nu <- prior.h$nu
    Psi <- prior.h$Psi
    a <- prior.h$a
    b <- prior.h$b
    aalpha <- prior.h$aalpha
    balpha <- prior.h$balpha

    if(is.na(L)) {
        L <- 10 
    } else { 
        if(length(L) != 1) {
            stop(paste0("L must be a constant"))
        }
    }
    if(all(is.na(m0))) {
        if(standardise) m0 <- rep(0, k)
        else m0 <- coefs.h
    } else { 
        if(length(m0) != k) {
            stop(paste0("'m0' must be a vector of length ", k))
        }
    }
    
    if(all(is.na(S0))){
        if(standardise) S0 <- 10*diag(k)
        else S0 <- cov.h
    } else { 
        if(!is.matrix(S0) | !all(dim(S0) == c(k, k))) {
            stop(paste0("'S0' must be a matrix of dimension ", k, "x", k))
        }
    }
    
    if(all(is.na(Psi))){
        if(standardise) Psi <- diag(k)
        else Psi <- 30*cov.h
    } else { 
        if(!is.matrix(Psi) | !all(dim(Psi) == c(k, k))) {
            stop(paste0("'Psi' must be a matrix of dimension ", k, "x", k))
        }
    }
    
    if(is.na(nu)) {
        nu <- k + 2
    } else{ 
        if(nu < k + 2){
            stop(paste0("'nu' must be larger than ", k + 2))
        }
    }
    
    if(is.na(a)) {
        a <- 2
    } else { 
        if(length(a) != 1) {
            stop(paste0("'a' must be a constant"))
        }
    }
    
    if(is.na(b)){
        if(standardise) b <- 2
        else b <- var.h
    } else { 
        if(length(b) != 1) {
            stop(paste0("'b' must be a constant"))
        }
    }
    
    if(L > 1) {
        if(is.na(aalpha)) {
            aalpha <- 2 
        } else { 
            if(length(aalpha) != 1) {
                stop(paste0("aalpha must be a constant"))
            }
        }
        
        if(is.na(balpha)) {
            balpha <- 2 
        } else {
            if(length(balpha) != 1) {
                stop(paste0("balpha must be a constant"))
            }
        }
    }
    
    if(L > 1){
        res0 <- bddp(y = data.h.marker,
        X = X0,
        prior = list(m0 = m0,
        S0 = S0,
        nu = nu,
        Psi = Psi,
        a = a,
        b = b,
        aalpha = aalpha,
        balpha = balpha,
        L = L),
        mcmc = mcmc,
        standardise = standardise)
    }
    
    if(L == 1){
        res0 <- regnth(y = data.h.marker,
        X = X0,
        prior = list(m0 = m0,
        S0 = S0,
        nu = nu,
        Psi = Psi,
        a = a,
        b = b),
        mcmc = mcmc,
        standardise = standardise)
    }
    if(density$compute) {
        if(!all(is.na(density$newdata)) && !inherits(density$newdata, "data.frame"))
            stop("Newdata (argument density) must be a data frame")
        if(!all(is.na(density$newdata)) && length(names.cov) != 0 && sum(is.na(match(names.cov, names(density$newdata)))))
            stop("Not all needed variables are supplied in newdata (argument density)")
        if(all(is.na(density$newdata))) {
            newdata <- cROCData(data.new, names.cov, group)
        } else {
            newdata <- na.omit(density$newdata[,names.cov,drop=FALSE])
        }
        density$Xhp <- predict(MM0, newdata = newdata)$X
        npred <- nrow(newdata)
        
        if(all(is.na(density$grid.h))) {
            density$grid.h <- seq(min(data.h.marker) - 1, max(data.h.marker) + 1, len = 200)
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
                    parallel::mclapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, L = L, yd = data.d.marker, Xd = X1, pauc = pauc, density = density, p = p, mc.cores = ncpus)
                } else if (do_snow) {                
                    if (is.null(cl)) {
                        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                        if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                            parallel::clusterSetRNGStream(cl)
                        }
                        res <- parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, L = L, yd = data.d.marker, Xd = X1, pauc = pauc, density = density, p = p)
                        parallel::stopCluster(cl)
                        res
                    } else {
                        if(!inherits(cl, "cluster")) {
                            stop("Class of object 'cl' is not correct")
                        } else {
                            parallel::parLapply(cl, seq_len(mcmc$nsave), doMCMCROC, res0 = res0, L = L, yd = data.d.marker, Xd = X1, pauc = pauc, density = density, p = p)
                        }                        
                    }
                }
            } else {
                lapply(seq_len(mcmc$nsave), doMCMCROC, res0 = res0, L = L, yd = data.d.marker, Xd = X1, pauc = pauc, density = density, p = p)
            }

        resBoot <- simplify2array(resBoot)    
        arocbbddp <- simplify2array(resBoot["ROC",])
        aucddp <- unlist(resBoot["AUC",])
        if(pauc$compute){
            paucddp <- unlist(resBoot["pAUC",])
        }
        if(density$compute){
            dens.h <- aperm(simplify2array(resBoot["dens.h",]), c(1,3,2))
            meanfun.h <- simplify2array(resBoot["meanfun.h",])

            meanfun.h.m <- apply(meanfun.h, 1, mean)
            meanfun.h.l <- apply(meanfun.h, 1, quantile, prob = 0.025)
            meanfun.h.h <- apply(meanfun.h, 1, quantile, prob = 0.975)
        }

    } else {
        stop("nsave should be larger than zero.")
    }

    AROC <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
    AROC[,1] <- apply(arocbbddp, 1,mean)
    AROC[,2] <- apply(arocbbddp, 1, quantile, prob = 0.025)
    AROC[,3] <- apply(arocbbddp, 1, quantile, prob = 0.975)
    AUC <- c(mean(aucddp), quantile(aucddp,c(0.025,0.975)))
    names(AUC) <- c("est","ql", "qh")

    res <- list()
    res$call <- match.call()
    res$data <- data
    res$missing.ind <- list(h = omit.h, d = omit.d)
    res$marker <- marker
    res$group <- group
    res$tag.h <- tag.h
    res$p <- p
    res$mcmc <- mcmc
    if(L > 1){
        res$prior <- list(m0 = m0,
        S0 = S0,
        nu = nu,
        Psi = Psi,
        a = a,
        b = b,
        aalpha = aalpha,
        balpha = balpha,
        L = L)
    }
    if(L == 1){
        res$prior <- list(m0 = m0,
        S0 = S0,
        nu = nu,
        Psi = Psi,
        a = a,
        b = b,
        L = L)
    }
    res$ROC <- AROC
    res$AUC <- AUC
    if(pauc$compute) {
        aux <- c(mean(paucddp), quantile(paucddp, c(0.025, 0.975)))
        if(pauc$focus == "FPF"){
            res$pAUC <- aux/pauc$value
        } else{
            res$pAUC <- aux/(1 - pauc$value)
        }
        names(res$pAUC) <- c("est","ql", "qh")
        attr(res$pAUC, "value") <- pauc$value
        attr(res$pAUC, "focus") <- pauc$focus
    }
    if(density$compute) {
        res$newdata <- newdata
        res$reg.fun.h <- data.frame(est = meanfun.h.m, ql = meanfun.h.l, qh = meanfun.h.h)
        res$dens <- list(grid = density$grid.h, dens = dens.h)
    }
    if(compute.lpml | compute.WAIC | compute.DIC) {
        term <- inf_criteria(y = data.h.marker, X = X0, res = res0)
    }
    
    if(compute.lpml) {
        if(L > 1){
            res$lpml <- lpml(y = data.h.marker, X = X0, res = res0, L = L, termsum = term)
        }
        if(L == 1){
            res$lpml <- lpmlp(y = data.h.marker, X = X0, res = res0, term = term)
        }
    }
    
    if(compute.WAIC) {
        if(L > 1){
            res$WAIC <- waicnp(y = data.h.marker, X = X0, res = res0, L = L, termsum = term)
        }
        if(L == 1){
            res$WAIC <- waicp(y = data.h.marker, X = X0, res = res0, term = term)
        }
    }
    
    if(compute.DIC) {
        if(L > 1){
            res$DIC <- dic(y = data.h.marker, X = X0, res = res0, L = L, termsum = term)
        }
        if(L == 1){
            res$DIC <- dicp(y = data.h.marker, X = X0, res = res0, term = term)
        }
    }
    
    # Results of the fit in the healthy population (neeeded to calculate predictive checks or other statistics)
    if(L > 1){
        res$fit <- list(formula = formula.h, mm = MM0, beta = res0$Beta, sd = sqrt(res0$Sigma), probs = res0$P)
    }
    if(L == 1){
        res$fit <- list(formula = formula.h, mm = MM0, beta = res0$Beta, sd = sqrt(res0$Sigma))
    }
    res$data_model <- list(y = list(h = data.h.marker, d = data.d.marker), X = list(h = X0, d = X1))
    class(res) <- c("AROC.bnp", "AROC")
    res
}
