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
mcmc = mcmccontrol()) {
    
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
    names.cov <- all.vars(formula.h)[-1]
    
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
    
    # Construct design matrix in diseased population (based on healthy)
    X1 <- predict(MM0, data.d)$X
    k <-  ncol(X0)
    
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
    
    nsimf <- mcmc$nsave
    
    y1 <- data.d.marker
    
    if(L > 1){
        prob <- res0$P
        beta <- res0$Beta
        sd <- sqrt(res0$Sigma2)
    }
    if(L == 1){
        prob <- NULL
        beta <- res0$Beta
        sd <- sqrt(res0$Sigma2)
    }
    
    if(L > 1){
        udddp <- matrix(0, nrow = n1, ncol = nsimf)
        
        for(l in 1:nsimf) {
            #udddp[,l] <- 1 - apply(t(prob[l,]*t(pnorm(y1, mean = X1%*%t(beta[l,,]), sd = rep(sd[l,], each = length(y1))))),1, sum)
            udddp[,l] <- 1 - apply(t(prob[l,]*t(pnorm(y1, mean = tcrossprod(X1, beta[l,,]), sd = rep(sd[l,], each = length(y1))))),1, sum)
        }
        weights <- matrix(0, nrow = n1, ncol = nsimf)
        for(l in 1:nsimf) {
            aux1 <- rexp(n1,1)
            weights[,l] <- aux1/sum(aux1)
        }
        arocbbddp <- matrix(0, nrow = np, ncol = nsimf)
        aucddp <- numeric(nsimf)
        if(pauc$compute) {
            paucddp <- numeric(nsimf)
            if(pauc$focus == "FPF"){
                # Truncated pv
                tudddp <- matrix(pmin(pauc$value, udddp), nrow = n1)
            }
        }
        
        for(j in 1:np) {
            arocbbddp[j,] <- colSums(weights*(udddp<=p[j]))
        }
        
        aucddp <- 1 - colSums(weights*udddp)
        if(pauc$compute) {
            if(pauc$focus == "FPF"){
                paucddp <- pauc$value - colSums(weights*tudddp)
            } else{
                for(l in 1:nsimf){
                    arocbbddp[1,] <- 0; arocbbddp[np,] <- 1
                    rocapp <- approxfun(p, arocbbddp[,l], method = "linear")
                    p1 <- uniroot(function(x) {rocapp(x) - pauc$value}, interval = c(0, 1))$root
                    paucddp[l] <- integrate(rocapp, lower = p1, upper = 1,
                    stop.on.error = FALSE)$value - (1 - p1)*pauc$value
                }
            }
        }
    }
    
    if(L == 1){
        up <- matrix(0, nrow = n1, ncol = nsimf)
        for(k in 1:nsimf) {
            up[,k] = 1 - pnorm(data.d.marker, mean = X1%*%beta[k,], sd = sd[k])
        }
        
        weights <- matrix(0, nrow = n1, ncol = nsimf)
        for(l in 1:nsimf) {
            aux1 <- rexp(n1,1)
            weights[,l] <- aux1/sum(aux1)
        }
        
        arocbbddp <- matrix(0, nrow = np, ncol = nsimf)
        aucddp <- numeric(nsimf)
        
        for(j in 1:np) {
            arocbbddp[j,] <- colSums(weights*(up<=p[j]))
        }
        
        aucddp <- 1 - colSums(weights*up)
        
        if(pauc$compute) {
            paucddp <- numeric(nsimf)
            if(pauc$focus == "FPF"){
                # Truncated pv
                tup <- matrix(pmin(pauc$value, up), nrow = n1)
            }
        }
        
        if(pauc$compute) {
            if(pauc$focus == "FPF"){
                paucddp <- pauc$value - colSums(weights*tup)
            } else{
                for(l in 1:nsimf){
                    arocbbddp[1,] <- 0; arocbbddp[np,] <- 1
                    rocapp <- approxfun(p, arocbbddp[,l], method = "linear")
                    p1 <- uniroot(function(x) {rocapp(x) - pauc$value}, interval = c(0, 1))$root
                    paucddp[l] <- integrate(rocapp, lower = p1, upper = 1,
                    stop.on.error = FALSE)$value - (1 - p1)*pauc$value
                }
            }
        }
    }
    
    AROC <- matrix(0, ncol = 3, nrow = np, dimnames = list(1:np, c("est","ql", "qh")))
    AROC[,1] <- apply(arocbbddp, 1,mean)
    AROC[,2] <- apply(arocbbddp, 1, quantile, prob = 0.025)
    AROC[,3] <- apply(arocbbddp, 1, quantile, prob = 0.975)
    AUC <- c(mean(aucddp), quantile(aucddp,c(0.025,0.975)))
    names(AUC) <- c("est","ql", "qh")
    
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
        X0p <- predict(MM0, newdata = newdata)$X
        npred <- nrow(newdata)
        
        if(all(is.na(density$grid.h))) {
            grid.h <- seq(min(data.h.marker) - 1, max(data.h.marker) + 1, len = 200)
        } else {
            grid.h <- density$grid.h
        }
        
        dens.h <- array(0, c(length(grid.h), nsimf, npred))
        meanfun.h <- matrix(0, nrow = npred, ncol = nsimf)
        
        for(k in 1:nsimf){
            if(L == 1){
                mu.h <- X0p%*%beta[k,]
                meanfun.h[,k] = mu.h
                for(l in 1:npred){
                    dens.h[, k, l] <- dnorm(grid.h, mean = mu.h[l], sd = sd[k])
                }
            }
            if(L > 1){
                mu.h <- tcrossprod(X0p, beta[k,,]) #X0p%*%t(beta[k,,])
                for(l in 1:npred){
                    aux0 <- norMix(mu = c(mu.h[l,]), sigma = sd[k,], w = prob[k,])
                    dens.h[, k, l] <- dnorMix(grid.h, aux0)
                    meanfun.h[l,k] <- sum(prob[k,]*t(mu.h[l,]))
                }
            }
        }
        
        meanfun.h.m <- apply(meanfun.h, 1, mean)
        meanfun.h.l <- apply(meanfun.h, 1, quantile, prob = 0.025)
        meanfun.h.h <- apply(meanfun.h, 1, quantile, prob = 0.975)
    }
    
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
        res$dens <- list(grid = grid.h, dens = dens.h)
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
        res$fit <- list(formula = formula.h, mm = MM0, beta = beta, sd = sd, probs = prob)
    }
    if(L == 1){
        res$fit <- list(formula = formula.h, mm = MM0, beta = beta, sd = sd)
    }
    res$data_model <- list(y = list(h = data.h.marker, d = data.d.marker), X = list(h = X0, d = X1))
    class(res) <- c("AROC.bnp", "AROC")
    res
}
